//  SimplexLP - Project #13: Revised Simplex in natID
//  Student: Šehzada Sijarić (19964)
//
//  MPS.h
//  Writer + reader for the (fixed/free) MPS format - the standard input
//  format of COIN-OR CLP and virtually every LP solver. The benchmark
//  exports every synthetic instance as .mps so the *identical* problem can
//  be fed to CLP for the external comparison required by the proposal:
//
//      clp instance.mps -primalsimplex   (or -barrier)
//
//  Supported subset:
//      ROWS    N (objective), E, L, G
//      COLUMNS numeric entries (no MARKER/integrality)
//      RHS     one RHS set
//      BOUNDS  none needed: all variables default to  0 <= x  (our model)
//      RANGES  not supported (rejected on read)
//
//  The reader handles everything the writer emits plus L/G rows, which it
//  converts to standard form by appending slack/surplus columns - so Netlib
//  style inequality models within this subset load as well.
#pragma once

#include "LPProblem.h"
#include <fstream>
#include <sstream>
#include <map>
#include <iomanip>

namespace lp
{

// ---------------------------------------------------------------- writer --
inline bool writeMPS(const StandardLP& P, const std::string& path,
                     std::string* errMsg = nullptr)
{
    std::ofstream f(path);
    if (!f.is_open())
    {
        if (errMsg) *errMsg = "cannot open file for writing: " + path;
        return false;
    }

    f << "NAME          " << P.name << "\n";
    f << "ROWS\n";
    f << " N  COST\n";
    for (int i = 0; i < P.rows(); ++i)
        f << " E  " << P.rowNames[i] << "\n";

    f << "COLUMNS\n";
    f << std::setprecision(15);
    for (int j = 0; j < P.colsTotal(); ++j)
    {
        const std::string& cn = P.colNames[j];
        if (P.c[j] != 0.0)
            f << "    " << std::left << std::setw(10) << cn
              << std::setw(10) << "COST" << " " << P.c[j] << "\n";
        for (const SparseEntry& e : P.A.column(j))
            f << "    " << std::left << std::setw(10) << cn
              << std::setw(10) << P.rowNames[e.row] << " " << e.val << "\n";
    }

    f << "RHS\n";
    for (int i = 0; i < P.rows(); ++i)
        if (P.b[i] != 0.0)
            f << "    " << std::left << std::setw(10) << "RHS"
              << std::setw(10) << P.rowNames[i] << " " << P.b[i] << "\n";

    f << "ENDATA\n";
    return true;
}

// ---------------------------------------------------------------- reader --
inline bool readMPS(const std::string& path, StandardLP& out,
                    std::string* errMsg = nullptr)
{
    std::ifstream f(path);
    if (!f.is_open())
    {
        if (errMsg) *errMsg = "cannot open file: " + path;
        return false;
    }

    enum class Section { None, Rows, Columns, Rhs, Ranges, Bounds };
    Section sec = Section::None;

    std::string objName;
    std::vector<char>        rowType;   // 'E','L','G'
    std::vector<std::string> rowNames;
    std::map<std::string, int> rowIndex;

    std::vector<std::string>   colNames;
    std::map<std::string, int> colIndex;
    std::vector<double>        objCoef;
    std::vector<std::vector<std::pair<int,double>>> colEntries;
    std::map<int, double> rhs;

    std::string name = "MPSLP";
    std::string line;

    auto fail = [&](const std::string& msg)
    {
        if (errMsg) *errMsg = msg;
        return false;
    };

    while (std::getline(f, line))
    {
        if (line.empty() || line[0] == '*')
            continue;

        std::istringstream is(line);

        if (line[0] != ' ' && line[0] != '\t')
        {
            std::string head;
            is >> head;
            if      (head == "NAME")    { is >> name; sec = Section::None; }
            else if (head == "ROWS")    sec = Section::Rows;
            else if (head == "COLUMNS") sec = Section::Columns;
            else if (head == "RHS")     sec = Section::Rhs;
            else if (head == "RANGES")  return fail("RANGES section is not supported");
            else if (head == "BOUNDS")  sec = Section::Bounds;
            else if (head == "ENDATA")  break;
            else return fail("unknown MPS section: " + head);
            continue;
        }

        switch (sec)
        {
            case Section::Rows:
            {
                std::string t, rn;
                is >> t >> rn;
                if (t == "N")
                {
                    if (objName.empty()) objName = rn; // first N row = objective
                }
                else if (t == "E" || t == "L" || t == "G")
                {
                    rowIndex[rn] = (int)rowNames.size();
                    rowNames.push_back(rn);
                    rowType.push_back(t[0]);
                }
                else
                    return fail("unsupported row type: " + t);
                break;
            }
            case Section::Columns:
            {
                std::string cn;
                is >> cn;
                if (cn == "MARKER" || line.find("'MARKER'") != std::string::npos)
                    return fail("integer MARKER sections are not supported");

                int j;
                auto it = colIndex.find(cn);
                if (it == colIndex.end())
                {
                    j = (int)colNames.size();
                    colIndex[cn] = j;
                    colNames.push_back(cn);
                    objCoef.push_back(0.0);
                    colEntries.emplace_back();
                }
                else
                    j = it->second;

                // one or two (row, value) pairs per line
                std::string rn; double v;
                while (is >> rn >> v)
                {
                    if (rn == objName)
                        objCoef[j] = v;
                    else
                    {
                        auto rit = rowIndex.find(rn);
                        if (rit == rowIndex.end())
                            return fail("COLUMNS references unknown row: " + rn);
                        colEntries[j].push_back({rit->second, v});
                    }
                }
                break;
            }
            case Section::Rhs:
            {
                std::string setName, rn; double v;
                is >> setName;
                while (is >> rn >> v)
                {
                    if (rn == objName)
                        continue; // objective constant - ignored
                    auto rit = rowIndex.find(rn);
                    if (rit == rowIndex.end())
                        return fail("RHS references unknown row: " + rn);
                    rhs[rit->second] = v;
                }
                break;
            }
            case Section::Bounds:
            {
                // only the trivial "LO ... 0" is acceptable in our subset
                std::string bt, setName, cn; double v = 0.0;
                is >> bt >> setName >> cn;
                is >> v;
                if (bt == "LO" && v == 0.0)
                    break;
                return fail("BOUNDS section beyond 'LO 0' is not supported");
            }
            default:
                break;
        }
    }

    if (rowNames.empty() || colNames.empty())
        return fail("MPS file has no rows or no columns");

    // ---- assemble standard form -------------------------------------------
    const int m  = (int)rowNames.size();
    const int n0 = (int)colNames.size();

    int nSlacks = 0;
    for (char t : rowType)
        if (t != 'E')
            ++nSlacks;

    out = StandardLP();
    out.name        = name;
    out.nStructural = n0;
    out.A.resize(m, n0 + nSlacks);
    out.b.assign(m, 0.0);
    out.c.assign(n0 + nSlacks, 0.0);
    out.colNames.resize(n0 + nSlacks);
    out.rowNames = rowNames;
    out.slackOfRow.assign(m, -1);

    for (int j = 0; j < n0; ++j)
    {
        out.colNames[j] = colNames[j];
        out.c[j]        = objCoef[j];
        for (auto& e : colEntries[j])
            out.A.add(e.first, j, e.second);
    }
    for (auto& kv : rhs)
        out.b[kv.first] = kv.second;

    int extra = n0;
    for (int i = 0; i < m; ++i)
    {
        if (rowType[i] == 'L')
        {
            out.A.add(i, extra, 1.0);
            out.colNames[extra] = "SLK" + std::to_string(i + 1);
            out.slackOfRow[i] = extra;
            ++extra;
        }
        else if (rowType[i] == 'G')
        {
            out.A.add(i, extra, -1.0);
            out.colNames[extra] = "SLK" + std::to_string(i + 1);
            ++extra;
        }
    }
    return true;
}

} // namespace lp
