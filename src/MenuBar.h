//  SimplexLP - Project #13: Revised Simplex in natID
//  Student: Šehzada Sijarić (19964)
//
//  MenuBar.h - three menus following the SDK example conventions:
//      App  : Settings | Quit
//      Run  : Solve | Step | Animate | Reset
//      Help : About
#pragma once
#include <gui/MenuBar.h>
#include "Constants.h"

class MenuBar : public gui::MenuBar
{
private:
    gui::SubMenu _subApp;
    gui::SubMenu _subRun;
    gui::SubMenu _subHelp;

protected:
    void populateAppMenu()
    {
        auto& items = _subApp.getItems();
        items[0].initAsActionItem(tr("settings"), cActSettings);
        items[1].initAsSeparator();
        items[2].initAsQuitAppActionItem(tr("Quit"), "q");
    }

    void populateRunMenu()
    {
        auto& items = _subRun.getItems();
        items[0].initAsActionItem(tr("solve"),   cActSolve,   "s");
        items[1].initAsActionItem(tr("step"),    cActStep,    "t");
        items[2].initAsActionItem(tr("animate"), cActAnimate, "p");
        items[3].initAsActionItem(tr("reset"),   cActReset,   "r");
    }

    void populateHelpMenu()
    {
        auto& items = _subHelp.getItems();
        items[0].initAsActionItem(tr("about"), cActAbout);
    }

public:
    MenuBar()
    : gui::MenuBar(3)
    , _subApp(cMenuApp,  tr("App"),  3)
    , _subRun(cMenuRun,  tr("Run"),  4)
    , _subHelp(cMenuHelp, tr("Help"), 1)
    {
        populateAppMenu();
        populateRunMenu();
        populateHelpMenu();
        _menus[0] = &_subApp;
        _menus[1] = &_subRun;
        _menus[2] = &_subHelp;
    }
};
