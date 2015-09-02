/*! \mainpage Borexino Montecarlo code
*
*\section sec_preamble Preamble
*
* The official Borexino Montecarlo simulation code is g4bx. To perform any reliable and complete simulation you must use g4bx.
*
* g4bx2 is under development. This means it is not tuned, debugged and complete so far.
 *
 *\section sec_intro_ Introduction
 *
*g4bx is the name of the Borexino Montecarlo code. 
*In 2014, a major revision of the code has been carried out. Consequently, a new version of the code (g4bx2) has been developed. This documentation shows details about g4bx2, how to install it and to use it. Many features are similar to the previous version of the code, so you can refer also to the guide <i>G4Bx Man Pages</i> made by D. Franco and I. Machulin in 2012.
*
* \section docs Documentation
*   Here you can find the links for some useful manuals:
*   - g4bx, <a href="https://bxweb.lngs.infn.it/docs/WorkingGroups/MC/g4bxman.pdf">G4Bx user manual</a> 
*   - bx_elec, <a href="https://bxweb.lngs.infn.it/docs/WorkingGroups/MC/UserManual_bxelec.pdf">bx_elec user manual</a> 
*   - Echidna, <a href="https://bxweb.lngs.infn.it/docs/WorkingGroups/MC/echidna_doc.pdf">Echidna user manual</a> 
*   - pdf version of this doxygen, <a href="https://bxweb.lngs.infn.it/docs/WorkingGroups/MC/doxycode.pdf">doxycode.pdf</a>
*
 * \section sec_start_ Getting ready
*\subsection get_the_code Getting the code
 *You can get the code from the usual cvs repository (<i>CVSROOT=username\@ui-borexino:/opt/exp_software/borexino/cvs</i>) typing
 *
 *<i>cvs co offline</i>
 *
 *then you have to compile the usual and the revised versions of the montecarlo, bx_elec and Echidna. 
 *
 *\subsection compile_g4bx Compile g4bx
 *Enter the folder <i>offline/bxmc/g4bx/setup/cnaf</i> and type
 *
 *<i> ./setup.sh </i>
 *
 *Now go back to <i>offline/bxmc/g4bx</i> and type
 *
 *<i>source please_source_me</i>
 *
 *and then
 *
 *<i>make</i>
 *
 * You are now ready to use g4bx.
 *
 *\subsection compile_g4bx2 Compile g4bx2
 *
 *Enter the folder <i>offline/bxmc/g4bx2/</i> and source the file <i>please_source_me</i> typing 
 *
 *<i>source please_source_me</i>
 *
 *Enter the folder <i>offline/bxmc/g4bx2/build</i> and type
 *
 <i>cmake -DGeant4_DIR=/opt/exp_software/borexino/geant4/geant4.10.00.p02/lib64/Geant4-10.0.2/ ../.</i>
 *
 *and then type
 *
 *<i>make</i>
 *
 *g4bx2 is now ready to be used. Remember to source  <i>please_source_me</i> every time you login at cnaf before using the code (eventually you could consider of adding it to your custom bashrc).
 *
 * If you want to use the old version of the montecarlo after having executed the new version, please run the script <i>unset_g4bx2.sh</i> in the folder <i>offline/bxmc/g4bx2/tools</i>. This is useful only if you run g4bx and g4bx2 by hand. If you use the automated script, this is not necessary.
 *
 *\subsection compile_bx_elec Compile bx_elec
 *
 *Enter the folder <i>offline/bxmc/bx_elec/setup/cnaf</i> and type
 *<i>./setup.sh</i>
 *
 *then move to the folder <i>offline/bxmc/bx_elec</i> and type
 *
 *<i>make</i>
 *
 *bx_elec is now ready to be used.
 *
 *\subsection compile_echidna Compile Echidna
 *
 *Enter the folder <i>offline/Echidna/setup/cnaf</i> and type
 *
 *<i>./setup.sh</i>
 *
 *source the 32-bit version of ROOT suggested by the output, by typing
 *
 *<i>. /opt/exp_software/borexino/root32/v5-34-24/bin/thisroot.sh</i>
 *
 *Now go back to <i>offline/Echidna</i> and type
 *
 *<i>make</i>
 *
 *Echidna is now ready to be used.
 *
 * \section sec_run_ Running g4bx2
 * \subsection sub_sec_vis Graphic mode
 * You can look at the geometry currently simulated by running g4bx2 without arguments. In fact, if in the build folder you just type
 *
 *<i>./g4bx2</i>
 *
 *a graphical window will open. The Qt graphical libraries let you interact with the visualization easily with the mouse. Just right-click on the
*picture and select "Mouse Actions" to decide whether you want to move, zoom or rotate the view.
*However, a few written commands could be useful. They have to be written in the "Session" field just below the picture. For example:
*
*<i>/vis/viewer/addCutawayPlane</i>
*
* cuts the detector with an YZ plane; 
*
*<i>/vis/viewer/set/viewpointVector -1 0 0</i>
*
*sets the pointview according to the vector (-1,0,0);
*
*<i>/vis/viewer/zoom 2</i>
*
*zooms the view of a factor of two.
*Many others can be found in the geant4 official guide.
*
*\subsection sub_sec_macro Macro file input
* You can execute the simulation by means of macro file (mac-file).
* They contain all the informations needed to carry on the simulation. To execute the simulation with a specific mac-file (e.g. <i>mac-file.mac</i>), simply type
*
*<i>./g4bx2 mac-file.mac</i>
*
* in the build folder.
*
*Some useful macro-files are already present in the build folder. For a detailed description of most of the mac-file options please refer to G4Bx man pages (they were not changed in the review of the code).
*
* \section sim_sub_ Automated script
*
*An automated script to send jobs in batch at cnaf is provided. Its name is <i>sim_sub_cnaf2.pl</i>. To use it, you have to copy this file from <i>offline/bxmc/g4bx2/tools</i> to <i>offline/.</i>. Be careful: in order to use successfully this automated script, you need g4bx2, bx_elec and Echidna to be correctly compiled and working. The first time you compile all the previous programs, you need to completely exit cnaf and enter again to be able to use the sim_sub script.
*
*To get the usage and a list of options, please type
*<i>./sim_sub_cnaf2.pl</i>
*
*Please note that, in order to execute the simulation with the revised montecarlo code, you have to set the option
*
*<i>version=2</i>
*
*\section g4rooter
*
* You have the possibility to create a rootfile directly as an output of Geant, instead of going through the chain of the electronics simulation and Echidna production.
* You can use the tool "g4rooter" located in <i>offline/bxmc/g4bx/tools</i> for g4bx and in <i>offline/bxmc/g4bx2/tools</i> for g4bx2.
* The procedure to make it work is the same, both for g4bx and g4bx2. Enter the folder <i>tools</i> in the <i>g4bx</i> (<i>g4bx2</i>) directory.
* In order to compile the tool, simply type
*
* <i>source compile_rooter</i>
*
* Now copy the binary file produced by Geant (.fil extension, e.g. output.fil) in the <i>tools</i> folder and type
*
* <i>./g4rooter output.fil</i>
*
* The result will be the root file output.root.
*
*\section mc_validation MC VALIDATION
*
* A tool created by S. Davini (stefano.davini@ge.infn.it) is provided to efficiently tune the MC on the basis of the calibration data.
* Up to now, this tool works only for the g4bx version of the code. To compile the tool, first of all link to the 32 bit version of root by typing
*
*<i>. /opt/exp_software/borexino/root32/v5-34-24/bin/thisroot.sh</i>
*
*Now enter the folder <i>offline/bxmc/g4bx/tools/mc_validation</i> and type
*
* <i>./compile_mc_validation.sh</i>
*
* Go back to the folder <i>offline</i> and run 
*
* <i>./mc_validation.pl</i>
*
* to display all the options you can provide to the script. First of all, you need to simulate the runs (simulate or fast_sim options). You can check if the simulation is over through the check option.
* When all the files are ready, you can run the mc_validation with the validation option. The output is the file mc_validation.root in the folder <i>offline/Echidna</i>.
* You can compare many outputs of the mc_validation through the script <i>compare_mc_validation</i> in the Echidna folder. In order to use it to compare mc_validation1.root, mc_validation2.root... type
*
* <i>./compare_mc_validation add:mc_validation1.root add:mc_validation2.root ...</i>
*
* \section contacts_ Contacts
* In case of problems contact:
*
* A. Caminata: alessio.caminata\@ge.infn.it
*
* S. Marcocci: simone.marcocci\@ge.infn.it
*/
