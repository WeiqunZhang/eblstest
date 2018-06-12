#include "MyTest.H"

#include <AMReX_MLEBABecLap.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_EB2.H>

using namespace amrex;

MyTest::MyTest ()
{
    readParameters();

    initGrids();

    initializeEB();

    initData();
}

void
MyTest::solve ()
{
#if 0
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

#endif

    for (int ilev = 0; ilev <= max_level; ++ilev) {
        amrex::VisMF::Write(phi[0], "phi-"+std::to_string(ilev));
        amrex::VisMF::Write(factory[ilev]->getVolFrac(), "vfrc-"+std::to_string(ilev));
    }
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
    acoef.resize(nlevels);
    bcoef.resize(nlevels);

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        dmap[ilev].define(grids[ilev]);
        const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
        const EB2::Level& eb_level = eb_is.getLevel(geom[ilev]);
        factory[ilev].reset(new EBFArrayBoxFactory(eb_level, geom[ilev], grids[ilev], dmap[ilev],
                                                   {2,2,2}, EBSupport::full));

        phi[ilev].define(grids[ilev], dmap[ilev], 1, 1, MFInfo(), *factory[ilev]);
        rhs[ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        acoef[ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            bcoef[ilev][idim].define(amrex::convert(grids[ilev],IntVect::TheDimensionVector(idim)),
                                     dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        }

        phi[ilev].setVal(0.0);
        rhs[ilev].setVal(0.0);
        acoef[ilev].setVal(0.0);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            bcoef[ilev][idim].setVal(1.0);
        }
    }
}

