//  SimplexLP - Project #13: Revised Simplex in natID
//  Student: Šehzada Sijarić (19964)
//
//  MainView.h
//  Hosts the three application tabs inside a gui::StandardTabView:
//
//      [ 2D visualization | 3D visualization | Benchmark ]
//
//  and routes the Run-menu / toolbar actions (Solve, Step, Animate, Reset)
//  to whichever visualization tab is currently active.
#pragma once

#include <gui/StandardTabView.h>
#include "viz/Tab2D.h"
#include "viz/Tab3D.h"
#include "bench/TabBench.h"

class MainView : public gui::StandardTabView
{
    Tab2D    _tab2D;
    Tab3D    _tab3D;
    TabBench _tabBench;

public:
    MainView()
    {
        addView(&_tab2D,    tr("tab2D"));
        addView(&_tab3D,    tr("tab3D"));
        addView(&_tabBench, tr("tabBench"));
        setCurrentViewPos(0);
    }

    // ---- routing of menu / toolbar actions ---------------------------------
    void solveActive()
    {
        switch (getCurrentViewPos())
        {
            case 0: _tab2D.solve(); break;
            case 1: _tab3D.solve(); break;
            default: break; // benchmark tab has its own Run button
        }
    }

    void stepActive()
    {
        switch (getCurrentViewPos())
        {
            case 0: _tab2D.step(); break;
            case 1: _tab3D.step(); break;
            default: break;
        }
    }

    void animateActive()
    {
        switch (getCurrentViewPos())
        {
            case 0: _tab2D.togglePlay(); break;
            case 1: _tab3D.togglePlay(); break;
            default: break;
        }
    }

    void resetActive()
    {
        switch (getCurrentViewPos())
        {
            case 0: _tab2D.reset(); break;
            case 1: _tab3D.reset(); break;
            default: break;
        }
    }
};
