/**
BSD-3-Clause
Copyright 2019 Alliance for Sustainable Energy, LLC
Redistribution and use in source and binary forms, with or without modification, are permitted provided
that the following conditions are met :
1.	Redistributions of source code must retain the above copyright notice, this list of conditions
and the following disclaimer.
2.	Redistributions in binary form must reproduce the above copyright notice, this list of conditions
and the following disclaimer in the documentation and/or other materials provided with the distribution.
3.	Neither the name of the copyright holder nor the names of its contributors may be used to endorse
or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT HOLDER, CONTRIBUTORS, UNITED STATES GOVERNMENT OR UNITED STATES
DEPARTMENT OF ENERGY, NOR ANY OF THEIR EMPLOYEES, BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
OR CONSEQUENTIAL DAMAGES(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <limits>
#include <cmath>
#include <string>

#include "emscripten.h"
#include "6par_jacobian.h"
#include "6par_lu.h"
#include "6par_search.h"
#include "6par_newton.h"
#include "6par_gamma.h"
#include "6par_solve.h"

extern "C" {
EMSCRIPTEN_KEEPALIVE
int sixparsolve(
    const char* celltype,
    double Vmp,
    double Imp,
    double Voc,
    double Isc,
    double bVoc, // beta_voc
    double aIsc, // alpha_isc
    double gPmp, // gamma_pmp,
    int Nser,
    double Tref,

    double &a,
    double &Il,
    double &Io,
    double &Rs,
    double &Rsh,
    double &Adj
) {
    int tech_id = module6par::monoSi;
    std::string stype = celltype;

    if (stype.find("mono") != std::string::npos) tech_id = module6par::monoSi;
    else if (stype.find("multi") != std::string::npos || stype.find("poly") != std::string::npos) tech_id = module6par::multiSi;
    else if (stype.find("cis") != std::string::npos) tech_id = module6par::CIS;
    else if (stype.find("cigs") != std::string::npos) tech_id = module6par::CIGS;
    else if (stype.find("cdte") != std::string::npos) tech_id = module6par::CdTe;
    else if (stype.find("amor") != std::string::npos) tech_id = module6par::Amorphous;
    else
        return 1; //could not determine cell type (mono,multi,cis,cigs,cdte,amorphous)

    module6par m( tech_id, Vmp, Imp, Voc, Isc, bVoc, aIsc, gPmp, Nser, Tref+273.15 );
    int err = m.solve_with_sanity_and_heuristics<double>(300,1e-7);
    if (err < 0)
        return 2; // could not solve, check inputs

    a = m.a;
    Il = m.Il;
    Io = m.Io;
    Rs = m.Rs;
    Rsh = m.Rsh;
    Adj = m.Adj;

    return 0;
}

}
