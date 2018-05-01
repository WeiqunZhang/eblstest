#ifndef MY_TEST_H_
#define MY_TEST_H_

#include <string>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_EBFabFactory.H>

class MyTest
{
public:
    MyTest ();
    ~MyTest ();

    void test ();

private:
    void readParameters ();
    void initGrids ();
    void initializeEBIS ();
    void initData ();

    int n_cell = 128;
    int max_grid_size = 64;
    std::string geom_type = "sphere";

    amrex::Geometry geom;
    amrex::BoxArray grids;
    amrex::DistributionMapping dmap;
    std::unique_ptr<amrex::EBFArrayBoxFactory> factory;

    amrex::Vector<std::string> moment_names;
    amrex::MultiFab moments;
};

#endif
