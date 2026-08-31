//  SimplexLP - Project #13: Revised Simplex in natID
//  Student: Šehzada Sijarić (19964)
//
//  View2DCanvas.h
//  gui::Canvas that renders a 2-variable LP:
//
//    * axes with ticks, constraint boundary lines (labelled C1, C2, ...)
//    * the feasible region polygon (Sutherland-Hodgman, see Polytope.h)
//    * the revised-simplex pivot path along the vertices, animated by a
//      cursor (current vertex highlighted, visited part drawn solid)
//    * the interior-point iterate path through the interior (dashed)
//    * the optimum marker with the objective value
//    * the objective gradient arrow and the iso-profit line through the
//      current vertex
//
//  Pan: primary-button drag.  Zoom: pinch / ctrl-wheel (onZoom).
//  All mapping world->screen is done manually so text and line widths stay
//  crisp at every zoom level.
#pragma once

#include <gui/Canvas.h>
#include <gui/Shape.h>
#include <gui/DrawableString.h>
#include "VizModel.h"

class View2DCanvas : public gui::Canvas
{
protected:
    lp::VizModel* _model = nullptr;

    gui::Size _size {600, 400};
    gui::Point _pan {0, 0};
    double _userZoom = 1.0;
    gui::Point _dragStart {0, 0};
    bool _showIPM = true;

public:
    explicit View2DCanvas(lp::VizModel* model)
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
        _pan = {0, 0};
        _userZoom = 1.0;
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
        _pan.x += p.x - _dragStart.x;
        _pan.y += p.y - _dragStart.y;
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

    // ---- world <-> screen ------------------------------------------------
    double scalePx() const
    {
        const double ext = (_model && _model->extent > 0) ? _model->extent : 10.0;
        const double avail = std::min(_size.width, _size.height) - 70.0;
        return std::max(10.0, avail) / ext * _userZoom;
    }

    gui::Point W(double wx, double wy) const
    {
        const double s = scalePx();
        return { 50.0 + _pan.x + wx * s,
                 _size.height - 45.0 + _pan.y - wy * s };
    }
    gui::Point W(const lp::Pt2& p) const { return W(p.x, p.y); }

    // ---- drawing ----------------------------------------------------------
    void onDraw(const gui::Rect& /*rect*/) override
    {
        if (!_model)
            return;

        drawAxes();

        if (!_model->hasProblem)
        {
            gui::DrawableString::draw(tr("hint2D"), gui::Point(60, 40),
                                      gui::Font::ID::SystemNormal, td::ColorID::Gray);
            return;
        }

        drawFeasibleRegion();
        drawConstraints();
        drawObjective();
        if (_showIPM)
            drawIPMPath();
        drawSimplexPath();
        drawOptimum();
    }

    void drawAxes()
    {
        const gui::Point o  = W(0, 0);
        const td::ColorID ax = td::ColorID::SysText;

        gui::Shape::drawLine({0, o.y}, {(double)_size.width, o.y}, ax, 1.2f);
        gui::Shape::drawLine({o.x, 0}, {o.x, (double)_size.height}, ax, 1.2f);

        // ticks every "nice" world unit
        const double ext  = _model ? _model->extent : 10.0;
        const double step = niceStep(ext / 5.0);
        for (double t = step; t <= ext * 4.0; t += step)
        {
            gui::Point px = W(t, 0);
            gui::Point py = W(0, t);
            if (px.x < _size.width - 8)
            {
                gui::Shape::drawLine({px.x, o.y - 3}, {px.x, o.y + 3}, ax, 1.0f);
                td::String lbl;
                lbl.format("%g", t);
                gui::DrawableString::draw(lbl, gui::Point(px.x - 8, o.y + 6),
                                          gui::Font::ID::SystemSmaller, td::ColorID::Gray);
            }
            if (py.y > 8)
            {
                gui::Shape::drawLine({o.x - 3, py.y}, {o.x + 3, py.y}, ax, 1.0f);
                td::String lbl;
                lbl.format("%g", t);
                gui::DrawableString::draw(lbl, gui::Point(o.x - 30, py.y - 7),
                                          gui::Font::ID::SystemSmaller, td::ColorID::Gray);
            }
        }
        gui::DrawableString::draw("x", gui::Point(_size.width - 16.0, o.y - 22),
                                  gui::Font::ID::SystemBold, td::ColorID::SysText);
        gui::DrawableString::draw("y", gui::Point(o.x + 8, 4),
                                  gui::Font::ID::SystemBold, td::ColorID::SysText);
    }

    void drawFeasibleRegion()
    {
        const auto& poly = _model->polygon;
        if (poly.size() < 3)
        {
            if (_model->hasProblem && poly.empty())
                gui::DrawableString::draw(tr("emptyRegion"), gui::Point(60, 60),
                                          gui::Font::ID::SystemBold, td::ColorID::Crimson);
            return;
        }

        std::vector<gui::Point> pts;
        pts.reserve(poly.size());
        for (const lp::Pt2& p : poly)
            pts.push_back(W(p));

        gui::Shape region;
        region.createPolygon(pts.data(), pts.size(), 1.6f);
        region.drawFillAndWire(td::ColorID::LightCyan, td::ColorID::SteelBlue);

        // region vertices
        for (const gui::Point& p : pts)
            drawMarker(p, 3.0, td::ColorID::SteelBlue);
    }

    void drawConstraints()
    {
        const double L = _model->extent * 4.0;
        int idx = 0;
        for (const auto& r : _model->user.rows)
        {
            ++idx;
            const double a0 = r.a.size() > 0 ? r.a[0] : 0.0;
            const double a1 = r.a.size() > 1 ? r.a[1] : 0.0;

            // two points of the line a0*x + a1*y = b inside [−L, L]
            lp::Pt2 p1, p2;
            if (std::fabs(a1) > std::fabs(a0))
            {
                p1 = {-0.1 * L, (r.b - a0 * (-0.1 * L)) / a1};
                p2 = { L,       (r.b - a0 * L) / a1};
            }
            else if (std::fabs(a0) > 1e-12)
            {
                p1 = {(r.b - a1 * (-0.1 * L)) / a0, -0.1 * L};
                p2 = {(r.b - a1 * L) / a0,           L};
            }
            else
                continue;

            const td::ColorID col = (r.rel == lp::Relation::Equal)
                                  ? td::ColorID::DarkOrange : td::ColorID::Gray;
            gui::Shape::drawLine(W(p1), W(p2), col, 1.0f, td::LinePattern::Dash);

            // label at the visible midpoint
            lp::Pt2 mid {(p1.x + p2.x) * 0.35, (p1.y + p2.y) * 0.35};
            td::String lbl;
            lbl.format("C%d", idx);
            gui::DrawableString::draw(lbl, W(mid), gui::Font::ID::SystemSmaller, col);
        }
    }

    void drawObjective()
    {
        // improvement direction in user space: +c for max, -c for min
        const auto& c = _model->user.c;
        if (c.size() < 2)
            return;
        double gx = c[0], gy = c[1];
        if (!_model->user.maximize)
        {
            gx = -gx;
            gy = -gy;
        }
        const double len = std::hypot(gx, gy);
        if (len < 1e-12)
            return;
        gx /= len;
        gy /= len;

        // arrow from the current vertex (or region centroid)
        lp::Pt2 base {_model->extent * 0.15, _model->extent * 0.15};
        if (const lp::IterSnapshot* s = _model->currentSnap())
            base = { s->x.size() > 0 ? s->x[0] : 0.0,
                     s->x.size() > 1 ? s->x[1] : 0.0 };

        const double aLen = _model->extent * 0.22;
        const lp::Pt2 tip {base.x + gx * aLen, base.y + gy * aLen};
        drawArrow(W(base), W(tip), td::ColorID::ForestGreen);

        // iso-profit line through the current vertex: c^T x = const
        const lp::Pt2 d {-gy, gx}; // perpendicular
        const double iso = _model->extent * 0.9;
        gui::Shape::drawLine(W(base.x - d.x * iso, base.y - d.y * iso),
                             W(base.x + d.x * iso, base.y + d.y * iso),
                             td::ColorID::ForestGreen, 1.0f, td::LinePattern::DashDot);
    }

    void drawSimplexPath()
    {
        if (!_model->hasPath())
            return;

        const int upto = _model->cursor;
        for (int k = 1; k < _model->pathLen(); ++k)
        {
            const lp::Pt2 a = _model->snap2(k - 1);
            const lp::Pt2 b = _model->snap2(k);
            const bool visited   = (k <= upto);
            const bool isPhase1  = (_model->simplex.path[(size_t)k].phase == 1);
            const td::ColorID col = isPhase1 ? td::ColorID::DarkOrange : td::ColorID::Crimson;

            if (visited)
                gui::Shape::drawLine(W(a), W(b), col, 3.0f);
            else
                gui::Shape::drawLine(W(a), W(b), td::ColorID::LightGray, 1.4f,
                                     td::LinePattern::Dash);
        }

        for (int k = 0; k <= upto && k < _model->pathLen(); ++k)
            drawMarker(W(_model->snap2(k)), 4.0, td::ColorID::Crimson);

        // highlight the current vertex
        if (const lp::IterSnapshot* s = _model->currentSnap())
        {
            const gui::Point p = W(_model->snap2(_model->cursor));
            drawRing(p, 8.0, td::ColorID::Crimson);
            td::String lbl;
            lbl.format("it %d  obj=%.4g", s->iter,
                       _model->std_.userObjective(s->x.data()));
            gui::DrawableString::draw(lbl, gui::Point(p.x + 10, p.y - 18),
                                      gui::Font::ID::SystemSmaller, td::ColorID::Crimson);
        }
    }

    void drawIPMPath()
    {
        const auto& path = _model->ipm.path;
        if (path.size() < 2)
            return;
        for (size_t k = 1; k < path.size(); ++k)
        {
            gui::Shape::drawLine(W(_model->ipmSnap2((int)k - 1)),
                                 W(_model->ipmSnap2((int)k)),
                                 td::ColorID::DarkMagenta, 1.6f, td::LinePattern::Dot);
        }
        for (size_t k = 0; k < path.size(); ++k)
            drawMarker(W(_model->ipmSnap2((int)k)), 2.5, td::ColorID::DarkMagenta);

        gui::DrawableString::draw(tr("lblIPMPath"),
                                  W(_model->ipmSnap2((int)path.size() / 2)),
                                  gui::Font::ID::SystemSmaller, td::ColorID::DarkMagenta);
    }

    void drawOptimum()
    {
        if (_model->simplex.status != lp::Status::Optimal)
            return;
        const lp::Pt2 opt { _model->simplex.x.size() > 0 ? _model->simplex.x[0] : 0.0,
                            _model->simplex.x.size() > 1 ? _model->simplex.x[1] : 0.0 };
        const gui::Point p = W(opt);
        drawMarker(p, 6.0, td::ColorID::Gold);
        drawRing(p, 10.0, td::ColorID::DarkOrange);

        td::String lbl;
        lbl.format("opt (%g, %g)  z=%.6g", opt.x, opt.y,
                   _model->std_.userObjective(&_model->simplex.x[0]));
        gui::DrawableString::draw(lbl, gui::Point(p.x + 12, p.y + 4),
                                  gui::Font::ID::SystemBold, td::ColorID::DarkOrange);
    }

    // ---- small helpers ------------------------------------------------------
    void drawMarker(const gui::Point& p, double r, td::ColorID col)
    {
        gui::Shape s;
        gui::Circle c(p, r);
        s.createCircle(c, 1.0f);
        s.drawFill(col);
    }

    void drawRing(const gui::Point& p, double r, td::ColorID col)
    {
        gui::Shape s;
        gui::Circle c(p, r);
        s.createCircle(c, 2.0f);
        s.drawWire(col);
    }

    void drawArrow(const gui::Point& a, const gui::Point& b, td::ColorID col)
    {
        gui::Shape::drawLine(a, b, col, 2.2f);
        const double dx = b.x - a.x, dy = b.y - a.y;
        const double L = std::hypot(dx, dy);
        if (L < 1e-9)
            return;
        const double ux = dx / L, uy = dy / L;
        const gui::Point h1 {b.x - 10 * ux + 5 * uy, b.y - 10 * uy - 5 * ux};
        const gui::Point h2 {b.x - 10 * ux - 5 * uy, b.y - 10 * uy + 5 * ux};
        gui::Shape::drawLine(b, h1, col, 2.2f);
        gui::Shape::drawLine(b, h2, col, 2.2f);
    }

    static double niceStep(double raw)
    {
        const double p = std::pow(10.0, std::floor(std::log10(std::max(raw, 1e-12))));
        const double m = raw / p;
        if (m < 1.5) return p;
        if (m < 3.5) return 2 * p;
        if (m < 7.5) return 5 * p;
        return 10 * p;
    }
};
