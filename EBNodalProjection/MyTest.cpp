#include "MyTest.H"

#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_EBIndexSpace.H>
#include <AMReX_EBTower.H>

using namespace amrex;

MyTest::MyTest ()
{
    readParameters();

    initGrids();

    initializeEBIS();

    EBTower::Build();
    AMReX_EBIS::reset();

    initData();
}

//
// Given vel, rhs & sig, this solves Div (sig * Grad phi) = Div vel + rhs.
// On return, vel becomes vel  - sig * Grad phi.
//
void
MyTest::solve ()
{
    static int count = 0;
    ++count;

    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;
    for (int idim = 0; idim < 2; ++idim) {
        if (Geometry::isPeriodic(idim)) {
            mlmg_lobc[idim] = LinOpBCType::Periodic;
            mlmg_hibc[idim] = LinOpBCType::Periodic;
        } else {
//            mlmg_lobc[idim] = LinOpBCType::Neumann;
//            mlmg_hibc[idim] = LinOpBCType::Neumann;
            mlmg_lobc[idim] = LinOpBCType::Dirichlet;
            mlmg_hibc[idim] = LinOpBCType::Dirichlet;
        }
    }
    mlmg_lobc[AMREX_SPACEDIM-1] = LinOpBCType::Neumann;
//    mlmg_lobc[AMREX_SPACEDIM-1] = LinOpBCType::Dirichlet;
    mlmg_hibc[AMREX_SPACEDIM-1] = LinOpBCType::Dirichlet;

    LPInfo info;
//    info.setAgglomeration(false);

    MLNodeLaplacian mlndlap(geom, grids, dmap, info, amrex::GetVecOfConstPtrs(factory));

    mlndlap.setHarmonicAverage(true);

    mlndlap.setGaussSeidel(false);

    if (sigma) {
        mlndlap.setCoarseningStrategy(MLNodeLaplacian::CoarseningStrategy::Sigma);
    }

    mlndlap.setDomainBC(mlmg_lobc, mlmg_hibc);

    for (int ilev = 0; ilev <= max_level; ++ilev) {
        mlndlap.setSigma(ilev, sig[ilev]);
    }

    for (int ilev = 0; ilev <= max_level; ++ilev) {
        amrex::VisMF::Write(vel[ilev], "velorig");
        amrex::VisMF::Write(factory[ilev]->getVolFrac(), "vfrc");
    }

    mlndlap.compRHS(amrex::GetVecOfPtrs(rhs), amrex::GetVecOfPtrs(vel), {}, {});

    for (int ilev = 0; ilev <= max_level; ++ilev) {
        amrex::VisMF::Write(rhs[ilev], "rhsorig");
    }

    MLMG mlmg(mlndlap);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(bottom_verbose);
    mlmg.setMaxIter(max_iter);
    mlmg.setMaxFmgIter(max_fmg_iter);

    // xxxxx
//    mlmg.setFinalSmooth(16);
//    mlmg.setBottomSmooth(8);

    // xxxxx
//    mlmg.setBottomSolver(MLMG::BottomSolver::smoother);
//    mlmg.setMaxIter(100);


    {
//        rhs[0].setVal(1.0);
    }


    Real mlmg_err = mlmg.solve(amrex::GetVecOfPtrs(phi), amrex::GetVecOfConstPtrs(rhs),
                               1.e-11, 0.0);

    mlndlap.updateVelocity(amrex::GetVecOfPtrs(vel), amrex::GetVecOfConstPtrs(phi));

    for (int ilev = 0; ilev <= max_level; ++ilev) {
        amrex::VisMF::Write(phi[ilev], "phi");
        amrex::VisMF::Write(vel[ilev], "vel");
    }

    mlndlap.compRHS(amrex::GetVecOfPtrs(rhs), amrex::GetVecOfPtrs(vel), {}, {});

    for (int ilev = 0; ilev <= max_level; ++ilev) {
        amrex::VisMF::Write(rhs[ilev], "rhs");
        amrex::Print() << "rhs.norm0() = " << rhs[ilev].norm0() << "\n";
        amrex::Print() << "rhs.norm1()/npoints = " << rhs[ilev].norm1() / grids[0].d_numPts() << "\n";
    }

    amrex::WriteSingleLevelPlotfile("plot", vel[0], {AMREX_D_DECL("xvel","yvel","zvel")}, geom[0], 0.0, 0);
}

void
MyTest::readParameters ()
{
    ParmParse pp;
    pp.query("max_level", max_level);
    pp.query("n_cell", n_cell);
    pp.query("max_grid_size", max_grid_size);

    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);
    pp.query("max_iter", max_iter);
    pp.query("max_fmg_iter", max_fmg_iter);
#ifdef AMREX_USE_HYPRE
    pp.query("use_hypre", use_hypre);
#endif

    pp.query("sigma", sigma);

    pp.query("geom_type", geom_type);
}

void
MyTest::initGrids ()
{
    int nlevels = max_level + 1;
    geom.resize(nlevels);
    grids.resize(nlevels);

    RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
    std::array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
//    std::array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1,0,0)};
    Geometry::Setup(&rb, 0, is_periodic.data());
    Box domain0(IntVect{AMREX_D_DECL(0,0,0)}, IntVect{AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1)});
    Box domain = domain0;
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        geom[ilev].define(domain);
        domain.refine(ref_ratio);
    }

    domain = domain0;
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        grids[ilev].define(domain);
        grids[ilev].maxSize(max_grid_size);
        domain.grow(-n_cell/4);   // fine level cover the middle of the coarse domain
        domain.refine(ref_ratio); 
    }
}

void
MyTest::initData ()
{
    int nlevels = max_level + 1;
    dmap.resize(nlevels);
    factory.resize(nlevels);
    phi.resize(nlevels);
    rhs.resize(nlevels);
    vel.resize(nlevels);
    sig.resize(nlevels);

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        dmap[ilev].define(grids[ilev]);
        factory[ilev].reset(new EBFArrayBoxFactory(geom[ilev], grids[ilev], dmap[ilev],
                                                   {2,2,2}, EBSupport::full));

        phi[ilev].define(amrex::convert(grids[ilev],IntVect::TheNodeVector()),
                         dmap[ilev], 1, 1, MFInfo(), *factory[ilev]);
        rhs[ilev].define(amrex::convert(grids[ilev],IntVect::TheNodeVector()),
                         dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        sig[ilev].define(grids[ilev], dmap[ilev], 1, 1, MFInfo(), *factory[ilev]);
        vel[ilev].define(grids[ilev], dmap[ilev], AMREX_SPACEDIM, 1, MFInfo(), *factory[ilev]);

        phi[ilev].setVal(0.0);
        sig[ilev].setVal(1.0);
        vel[ilev].setVal(0.0);
        const Real* dx = geom[ilev].CellSize();
        const Real h = dx[0];
        for (MFIter mfi(vel[ilev]); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            FArrayBox& fab = vel[ilev][mfi];
            int icomp = 0; // vx
            int ncomp = 1;
            fab.ForEachIV(bx, icomp, ncomp, [=] (Real& v_x, const IntVect& iv) {
                    Real rx = (iv[0]+0.5)*h - 0.5;
                    Real ry = (iv[1]+0.5)*h - 0.5;
                    Real r = std::sqrt(rx*rx+ry*ry);
                    Real fac = std::exp(-(r*r/(0.16*0.16)));
                    v_x = v_x + 2.0*r*ry/r*fac;
                });
            icomp = 1; // vy
            fab.ForEachIV(bx, icomp, ncomp, [=] (Real& v_y, const IntVect& iv) {
//                    if (iv[1] > n_cell/2+5) v_y = 1.0;
                    Real rx = (iv[0]+0.5)*h - 0.5;
                    Real ry = (iv[1]+0.5)*h - 0.5;
                    Real r = std::sqrt(rx*rx+ry*ry);
                    Real fac = std::exp(-(r*r/(0.16*0.16)));
                    v_y = v_y - 2.0*r*rx/r*fac;
                });
        }

        // xxxxx
//        vel[ilev].setVal(0.0, 0, 1, 0);
//        vel[ilev].setVal(1.0, 1, 1, 0);
    }

}

