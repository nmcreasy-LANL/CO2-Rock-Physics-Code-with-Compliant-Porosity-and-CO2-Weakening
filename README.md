# CO2-Rock-Physics-Code-with-Compliant-Porosity-and-CO2-Weakening
LANL software release number: O4934

An NRAP funded code development for CO2 storage. This code converts reservoir simulations to seismic properties for underground CO2 storage. 

Instructions:
python rock_physics_co2.py testdir mod1_090y.csv baseline_y0.csv

Csv files need to have the following columns:
Basline file - 
element, x, y, z, material ID, permfactor, porosity, horizontal permeability, vertical permeability, Pore pressure (Pascal), T (C), and salt_fractions

Time step files (e.g., mod1_090y.csv):
element, x, y, z, pore pressure, temp, gas saturation, salt

These csv files come from TOUGH (reservoir simulator) outputs. 




This program is Open-Source under the BSD-3 License.

 

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

 

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 

Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 

Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
