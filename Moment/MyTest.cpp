#include "MyTest.H"
#include <AMReX_ParmParse.H>
#include <AMReX_EBTower.H>
#include <AMReX_EBIndexSpace.H>
#include <AMReX_EBLevelGrid.H>
#include <AMReX_VisMF.H>

using namespace amrex;

MyTest::MyTest ()
{
    readParameters();

    initGrids();

    initializeEBIS();

    EBTower::Build();
//    AMReX_EBIS::reset();
    
    initData();
}

MyTest::~MyTest ()
{}

void
MyTest::readParameters ()
{
    ParmParse pp;
    pp.query("n_cell", n_cell);
    pp.query("max_grid_size", max_grid_size);
    pp.query("geom_type", geom_type);
}

void
MyTest::initGrids ()
{
    RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
    std::array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
    Geometry::Setup(&rb, 0, is_periodic.data());
    Box domain0(IntVect{AMREX_D_DECL(0,0,0)}, IntVect{AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1)});

    geom.define(domain0);

    grids.define(domain0);
    grids.maxSize(max_grid_size);

    dmap.define(grids);
}

void
MyTest::initData ()
{
    factory.reset(new EBFArrayBoxFactory(geom, grids, dmap, {2,2,2}, EBSupport::full));

    int ncomp = IndMomSpaceDim::size();
    moment_names.resize(ncomp);
    for(int ivec = 0; ivec < ncomp; ivec++)
    {
        IvSpaceDim iv = IndMomSpaceDim::getIndex(ivec);
        std::string xint = std::string("x^") + std::to_string(iv[0]);
        std::string yint = std::string(" y^") + std::to_string(iv[1]);
#if BL_SPACEDIM==3
        std::string zint = std::string(" z^") + std::to_string(iv[2]);
#endif
        std::string integrand = D_TERM(xint, +yint, +zint);
        std::string name = std::string("\\int ") + integrand + std::string(" dV") ;
        moment_names[ivec] = name;
        amrex::Print() << "i = " << ivec << ": " << name << "\n";
    }

    moments.define(grids, dmap, ncomp, 0, MFInfo(), *factory);

    EBLevelGrid eblg(grids, dmap, geom.Domain(), 1);
    const EBISLayout& ebisl = eblg.getEBISL();
    for (MFIter mfi(moments); mfi.isValid(); ++mfi)
    {
        moments[mfi].setVal(0);

        const EBISBox& ebisBox = ebisl[mfi];
        Box               grid = mfi.fabbox();
        IntVectSet ivs(grid);
        VoFIterator vofit(ivs, ebisBox.getEBGraph());
        const Vector<VolIndex>& vofs = vofit.getVector();
        for(int ivof = 0 ; ivof < vofs.size(); ivof++)
        {
            IndMomSpaceDim volmom = ebisBox.getEBData().getVolumeMoments(vofs[ivof]);
            for(MomItSpaceDim momit; momit.ok(); ++momit)
            {
                int ivar = IndMomSpaceDim::indexOf(momit());
                moments[mfi](vofs[ivof].gridIndex(), ivar) = volmom[momit()];
            }
        }
    }
}

void
MyTest::test ()
{
    VisMF::Write(moments, "moments");
}
