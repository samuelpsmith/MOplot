# MOplot

MOplos-MOproblems is a simple python script to plot and annotate molecular orbital energy diagrams from a csv file. It is meant to offer some improvements over alternative scripts, which may not annotate points. 
    
Copyright (C) 2023  Samuel Smith

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>

Contents:

1. MOplots_MOproblems.py is a simple python script to plot and annotate molecular orbital diagrams from a csv file.
2. MTPP.csv is an example csv file generated manually from electronic structure code output. Users may choose to generate such a file manually or programmatically. 

Quickstart:
1. Open MOplots_MOproblems and point to your file.
2. Choose between plotting simple categorical compound vs energy level, or try to use swarmplot to dodge/jitter the points to show degenerate energy levels side by side. I think the former looks better because the latter is asymmetric, and so it is enabled by default.
3. Run the program. A simple loop annotates energy levels with labels from the csv. Degeneracies of up to 4 are suppported. To change how close energies recognized as degenerate levels are, you can change the parameter degen, which is in eV. 
