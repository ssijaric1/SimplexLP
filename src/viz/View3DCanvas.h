//  SimplexLP - Project #13: Revised Simplex in natID
//  Student: Šehzada Sijarić (19964)
//
//  View3DCanvas.h
//  Interactive 3D view of a three-variable LP: the feasible polytope is
//  rendered with a classic painter's algorithm (faces sorted back-to-front
//  by view depth), the simplex pivot path runs along the polytope's edges
//  and the interior-point trajectory cuts through its interior.
//
//  Camera model (orthographic, plenty for a teaching tool):
//      r = Rx(pitch) * Rz(yaw) * p
//      screen.x = r.x,  screen.y = r.z (up),  depth = r.y (away from viewer)
//  Dragging rotates (yaw/pitch), the trackpad/mouse-wheel zoom scales.
//
//  Everything is drawn with gui::Shape / gui::DrawableString on a plain
//  gui::Canvas - no external 3D engine, as discussed in the proposal.
#pragma once

#include <gui/Canvas.h>
#include <gui/Shape.h>
#include <gui/DrawableString.h>
#include "VizModel.h"
#include <algorithm>
#include <cmath>

class View3DCanvas : public gui::Canvas
{
protected:
    lp::VizModel* _model = nullptr;

    gui::Size  _size {600, 400};
    double     _yaw   = 0.65;   // radians
    double     _pitch = 0.42;
    double     _userZoom = 1.0;
    gui::Point _dragStart {0, 0};
    bool       _showIPM = true;

public:
    explicit View3DCanvas(lp::VizModel* model)
    : gui::Canvas({gui::InputDevice::Event::PrimaryClicks,
                   gui::InputDevice::Event::CursorDrag,
                   gui::InputDevice::Event::Zoom})
    , _model(model)
    {
        enableResizeEvent(true);
    }

    void showIPM(bool b) { _showIPM = b; reDraw(); }

    void resetView()
    {
        _yaw = 0.65; _pitch = 0.42; _userZoom = 1.0;
        reDraw();
    }

protected:
    void onResize(const gui::Size& newSize) override { _size = newSize; }

    void onPrimaryButtonPressed(const gui::InputDevice& dev) override
    {
        _dragStart = dev.getFramePoint();
    }

    void onCursorDragged(const gui::InputDevice& dev) override
    {
        const gui::Point p = dev.getFramePoint();
        _yaw   += (p.x - _dragStart.x) * 0.010;
        _pitch += (p.y - _dragStart.y) * 0.010;
        const double lim = 1.55; // keep away from gimbal flip
        if (_pitch >  lim) _pitch =  lim;
        if (_pitch < -lim) _pitch = -lim;
        _dragStart = p;
        reDraw();
    }

    bool onZoom(const gui::InputDevice& dev) override
    {
        _userZoom *= dev.getScale();
        if (_userZoom < 0.05) _userZoom = 0.05;
        if (_userZoom > 50.0) _userZoom = 50.0;
        reDraw();
        return true;
    }

    // ---- model space -> rotated camera space -> screen ---------------------
    struct R3 { double x, y, z; };

    R3 rotate(const lp::Pt3& p) const
    {
        const double cy = std::cos(_yaw),   sy = std::sin(_yaw);
        const double cp = std::cos(_pitch), sp = std::sin(_pitch);
        // yaw about the model z-axis
        const double x1 =  cy * p.x - sy * p.y;
        const double y1 =  sy * p.x + cy * p.y;
        const double z1 =  p.z;
        // pitch about the screen-horizontal axis
        const double y2 =  cp * y1 - sp * z1;
        const double z2 =  sp * y1 + cp * z1;
        return {x1, y2, z2};
    }

    double drawScale() const
    {
        const double ext = (_model && _model->extent > 1e-9) ? _model->extent : 10.0;
        const double s = 0.34 * std::min((double)_size.width, (double)_size.height) / ext;
        return s * _userZoom;
    }

    gui::Point project(const R3& r) const
    {
        const double s  = drawScale();
        const double cx = _size.width  * 0.5;
        const double cy = _size.height * 0.58;
        // model is mostly in the positive octant; shift its center to view center
        const double ext = (_model && _model->extent > 1e-9) ? _model->extent : 10.0;
        const R3 c = rotate({ext * 0.5, ext * 0.5, ext * 0.4});
        return gui::Point(cx + (r.x - c.x) * s, cy - (r.z - c.z) * s);
    }

    gui::Point project(const lp::Pt3& p) const { return project(rotate(p)); }

    // -------------------------------------------------------------- drawing
    void onDraw(const gui::Rect& /*rect*/) override
    {
        if (!_model || !_model->hasProblem)
        {
            gui::DrawableString::draw(tr("hint3D"), gui::Point(40, 40),
                                      gui::Font::ID::SystemNormal, td::ColorID::Gray);
            return;
        }

        drawAxes();

        const lp::Polytope3& P = _model->poly3;
        if (P.empty())
        {
            gui::DrawableString::draw(tr("emptyRegion"), gui::Point(60, 60),
                                      gui::Font::ID::SystemBold, td::ColorID::Crimson);
            drawPaths();
            return;
        }

        // painter's algorithm: sort faces far -> near (larger depth first)
        std::vector<std::pair<double, int>> order;
        order.reserve(P.faces.size());
        for (int f = 0; f < (int)P.faces.size(); ++f)
        {
            const R3 c = rotate(P.faces[(size_t)f].centroid);
            order.push_back({c.y, f});
        }
        std::sort(order.begin(), order.end(),
                  [](const auto& a, const auto& b) { return a.first > b.first; });

        for (const auto& [depth, f] : order)
        {
            (void)depth;
            drawFace(P, P.faces[(size_t)f]);
        }

        drawPaths();
        drawOptimum();
        drawCurrentLabel();
    }

    void drawFace(const lp::Polytope3& P, const lp::Face3& face)
    {
        if (face.verts.size() < 3)
            return;

        const lp::Plane3& pl = P.planes[(size_t)face.plane];
        const bool isBoxCap = (pl.sourceRow < 0) && (pl.rhs > 1e-12); // safety box
        const bool isAxis   = (pl.sourceRow < 0) && !isBoxCap;        // x=0, y=0, z=0

        std::vector<gui::Point> pts;
        pts.reserve(face.verts.size());
        for (int vi : face.verts)
            pts.push_back(project(P.verts[(size_t)vi]));

        gui::Shape shape;
        if (isBoxCap)
        {
            // artificial cap of an unbounded region: wire only, very light
            shape.createPolygon(pts.data(), pts.size(), 0.8f,
                                td::LinePattern::Dash);
            shape.drawWire(td::ColorID::LightGray);
            return;
        }

        shape.createPolygon(pts.data(), pts.size(), 1.4f);
        if (isAxis)
            shape.drawFillAndWire(td::ColorID::LightYellow, td::ColorID::Gray);
        else
            shape.drawFillAndWire(td::ColorID::LightCyan, td::ColorID::SteelBlue);

        // label user constraints at their centroid
        if (pl.sourceRow >= 0)
        {
            td::String lbl;
            lbl.format("C%d", pl.sourceRow + 1);
            gui::DrawableString::draw(lbl, project(face.centroid),
                                      gui::Font::ID::SystemSmaller,
                                      td::ColorID::SteelBlue);
        }
    }

    void drawAxes()
    {
        const double ext = _model->extent * 1.25;
        const gui::Point o = project(lp::Pt3{0, 0, 0});
        const gui::Point px = project(lp::Pt3{ext, 0, 0});
        const gui::Point py = project(lp::Pt3{0, ext, 0});
        const gui::Point pz = project(lp::Pt3{0, 0, ext});

        gui::Shape::drawLine(o, px, td::ColorID::SysText, 1.2f);
        gui::Shape::drawLine(o, py, td::ColorID::SysText, 1.2f);
        gui::Shape::drawLine(o, pz, td::ColorID::SysText, 1.2f);

        gui::DrawableString::draw("x", px, gui::Font::ID::SystemBold, td::ColorID::SysText);
        gui::DrawableString::draw("y", py, gui::Font::ID::SystemBold, td::ColorID::SysText);
        gui::DrawableString::draw("z", pz, gui::Font::ID::SystemBold, td::ColorID::SysText);
    }

    void drawMarker(const gui::Point& p, double r, td::ColorID color)
    {
        gui::Shape s;
        s.createCircle(gui::Circle(p, r), 1.0f);
        s.drawFill(color);
    }

    void drawRing(const gui::Point& p, double r, td::ColorID color, float w)
    {
        gui::Shape s;
        s.createCircle(gui::Circle(p, r), w);
        s.drawWire(color);
    }

    void drawPaths()
    {
        // ---- interior-point trajectory (through the interior) -------------
        if (_showIPM && !_model->ipm.path.empty())
        {
            const int n = (int)_model->ipm.path.size();
            std::vector<gui::Point> pts;
            pts.reserve((size_t)n);
            for (int k = 0; k < n; ++k)
                pts.push_back(project(_model->ipmSnap3(k)));

            gui::Shape line;
            line.createPolyLine(pts.data(), pts.size(), 1.6f, td::LinePattern::Dot);
            line.drawWire(td::ColorID::DarkMagenta);
            for (const gui::Point& p : pts)
                drawMarker(p, 2.0, td::ColorID::DarkMagenta);
        }

        // ---- simplex pivot path (along the edges) --------------------------
        if (!_model->simplex.path.empty())
        {
            const int last = _model->cursor;
            const int n    = (int)_model->simplex.path.size();

            // future part of the path: light dashed preview
            if (last + 1 < n)
            {
                std::vector<gui::Point> rest;
                for (int k = last; k < n; ++k)
                    rest.push_back(project(_model->snap3(k)));
                gui::Shape line;
                line.createPolyLine(rest.data(), rest.size(), 1.2f,
                                    td::LinePattern::Dash);
                line.drawWire(td::ColorID::LightGray);
            }

            // visited part: solid, phase-aware coloring
            for (int k = 1; k <= last && k < n; ++k)
            {
                const gui::Point a = project(_model->snap3(k - 1));
                const gui::Point b = project(_model->snap3(k));
                const bool ph1 = _model->simplex.path[(size_t)k].phase == 1;
                gui::Shape::drawLine(a, b,
                                     ph1 ? td::ColorID::DarkOrange
                                         : td::ColorID::Crimson, 3.0f);
            }
            for (int k = 0; k <= last && k < n; ++k)
                drawMarker(project(_model->snap3(k)), 3.0, td::ColorID::Crimson);

            // current vertex
            drawRing(project(_model->snap3(last)), 7.0, td::ColorID::Crimson, 2.0f);
        }
    }

    void drawOptimum()
    {
        if (_model->simplex.status != lp::Status::Optimal ||
            _model->simplex.x.size() < 3)
            return;
        const lp::Pt3 opt {_model->simplex.x[0], _model->simplex.x[1], _model->simplex.x[2]};
        const gui::Point p = project(opt);
        drawMarker(p, 4.5, td::ColorID::Gold);
        drawRing(p, 8.0, td::ColorID::DarkOrange, 2.0f);

        td::String lbl;
        lbl.format("  opt (%.4g, %.4g, %.4g)", opt.x, opt.y, opt.z);
        gui::DrawableString::draw(lbl, gui::Point(p.x + 8, p.y - 14),
                                  gui::Font::ID::SystemBold, td::ColorID::DarkOrange);
    }

    void drawCurrentLabel()
    {
        const lp::IterSnapshot* s = _model->currentSnap();
        if (!s)
            return;
        td::String lbl;
        lbl.format("iter %d   phase %d   obj = %.6g",
                   s->iter, s->phase,
                   _model->std_.objSign * s->objective);
        gui::DrawableString::draw(lbl, gui::Point(12, 10),
                                  gui::Font::ID::SystemNormal, td::ColorID::SysText);
    }
};
