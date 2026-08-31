//  SimplexLP - Project #13: Revised Simplex in natID
//  Student: Šehzada Sijarić (19964)
//
//  MainWindow.h - top-level window: menu bar, toolbar and the tabbed main
//  view. Run-menu actions are routed to the active visualization tab;
//  App->Settings opens the language dialog; Help->About shows project info.
#pragma once
#include <gui/Window.h>
#include "MenuBar.h"
#include "ToolBar.h"
#include "MainView.h"
#include "DialogSettings.h"
#include "Constants.h"

class MainWindow : public gui::Window
{
protected:
    MenuBar  _mainMenuBar;
    ToolBar  _toolBar;
    MainView _mainView;

protected:
    void onInitialAppearance() override
    {
        _mainView.setFocus();
    }

    bool shouldClose() override
    {
        return true;
    }

    bool onActionItem(gui::ActionItemDescriptor& aiDesc) override
    {
        auto [menuID, firstSubMenuID, lastSubMenuID, actionID] = aiDesc.getIDs();
        (void)firstSubMenuID;
        (void)lastSubMenuID;

        switch (menuID)
        {
            case cMenuApp:
            {
                if (actionID == cActSettings)
                {
                    auto pDlg = getAttachedWindow(cSettingsDlgID);
                    if (pDlg)
                        pDlg->setFocus();
                    else
                    {
                        DialogSettings* pSettingsDlg =
                            new DialogSettings(this, cSettingsDlgID);
                        pSettingsDlg->keepOnTopOfParent();
                        pSettingsDlg->setMainTB(&_toolBar);
                        pSettingsDlg->open();
                    }
                    return true;
                }
                break;
            }
            case cMenuRun:
            {
                switch (actionID)
                {
                    case cActSolve:   _mainView.solveActive();   return true;
                    case cActStep:    _mainView.stepActive();    return true;
                    case cActAnimate: _mainView.animateActive(); return true;
                    case cActReset:   _mainView.resetActive();   return true;
                    default: break;
                }
                break;
            }
            case cMenuHelp:
            {
                if (actionID == cActAbout)
                {
                    showAlert(tr("aboutTitle"), tr("aboutText"));
                    return true;
                }
                break;
            }
            default:
                break;
        }
        return false;
    }

public:
    MainWindow()
    : gui::Window(gui::Geometry(40, 40, 1240, 800))
    {
        setTitle(tr("appTitle"));
        _mainMenuBar.setAsMain(this);
        setToolBar(_toolBar);
        setCentralView(&_mainView);
    }
};
