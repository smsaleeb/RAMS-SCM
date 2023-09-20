# RAMS-SCM
RAMS Single Column Microphysics Code

-------------------------------------------

Here we are running the SCM by using the RAMS model single column input/output.

Run "zz.test.all.sc" to perform the tests for all the sample simulations in "AA.scm.in.out.data".
This will ingest in the inputs from "scm.in", run the physics, and generate "scm.out" files.
These scm.in/*.txt and scm.out/*.txt files are not included in the git repository since they 
make change with version and there are a lot of them. These are generated in the "rams_scm_check"
workspace on my machine by running RAMS for a single timestep and creating the before (scm.in)
files and after (scm.out) conditions for each RAMS prognostic variable. These files are copied
to the "rams_git_scm" directory for operation. The "IN" files are read into the SCM, then the 
SCM is run 1 timestep and outputs the "OUT" files. Then we compare the RAMS "scm.out" files with
the SCM "scm.out" files and check for differences. Differences are output to "zz.out.diffs.txt".

-------------------------------------------

Note that for testing INITIALIZATION, we run the "y.zero.sc" script to zero out the aerosol
variables in the "scm.in" directory that are copied in here. This so we can truly test the
aerosol initialization. Within the "scm.in.0" directories, RAMS has already initializaed
the aerosols, so we zero this out for just the aerosol variables and leave the others alone
so we can test SCM aerosol initialization.

-------------------------------------------

For custom testing like using "test.sedimentation.RAMSIN" run "test.sedimentation.sh" 
or create something similar. Set flag "initcustom" in the RAMSIN. Then you need to custom
set your variables in the file "mic_init_scm.f90" subroutine "init_custom".
