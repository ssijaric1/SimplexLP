//  SimplexLP - Project #13: Revised Simplex in natID
//  Student: Šehzada Sijarić (19964)
//
//  ToolBar.h - quick access to Settings and Solve. The images come from the
//  SDK's Common/Icons set (declared in res/main.xml), the same resources the
//  SDK examples ship with.
#pragma once
#include <gui/ToolBar.h>
#include <gui/Image.h>
#include "Constants.h"

class ToolBar : public gui::ToolBar
{
    gui::Image _imgSettings;
    gui::Image _imgSolve;

public:
    ToolBar()
    : gui::ToolBar("mainTB", 2)
    , _imgSettings(":settings")
    , _imgSolve(":start")
    {
        addItem(tr("settings"), &_imgSettings, tr("settingsTT"),
                cMenuApp, 0, 0, cActSettings);
        addItem(tr("solve"), &_imgSolve, tr("solveTT"),
                cMenuRun, 0, 0, cActSolve);
    }
};
