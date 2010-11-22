--------------------------------
readme
--------------------------------

To compare synthetics and data, the following steps are recommended:

1. Make sure that both synthetic and observed seismograms have the
   correct station/event and timing information.

2. Convolve synthetic seismograms with a source time function with the
   half duration speciﬁed in the CMTSOLUTION ﬁle, provided, as recommended,
   you used a zero half duration in the SEM simulations.

3. Resample both observed and synthetic seismograms to a common sampling rate.

4. Cut the records using the same window.

5. Remove the trend and mean from the records and taper them.

6. Remove the instrument response from the observed seismograms (recommended)
   or convolve the synthetic seismograms with the instrument response.

7. Make sure that you apply the same ﬁlters to both observed and synthetic seismograms.
   Preferably, avoid ﬁltering your records more than once.

8. Now, you are ready to compare your synthetic and observed seismograms.


Example:
--------

- data processing:

  > process_data.pl -m CMTSOLUTION -s 1.0 -l 0/4000 -i -f -t 40/500 -p -x bp DATA/1999.330*.BX?.SAC

    - adds CMTSOLUTION information
    - resamples ( 1 s sampling rate)
    - cuts ( between 0 and 4000 s)
    - transfers data to displacement and filters (between 40 and 500 s)
    - picks P & S arrivals
    - adds suffix 'bp'

- synthetics processing:

  > process_syn.pl -m CMTSOLUTION -h -a STATIONS -s 1.0 -l 0/4000 -f -t 40/500 -p -x bp SEM/*.BX?.sem

    - adds CMTSOLUTION information
    - convolves with source-time function
    - adds STATIONS information
    - resamples ( 1 s sampling rate)
    - cuts ( between 0 and 4000 s)
    - filters (between 40 and 500 s) compatible with transfer function
    - picks P & S arrivals
    - adds suffix 'bp'

to rotate the horizontal components:

  > rotate.pl -l 0 -L 4000 -d DATA/*.BXE.SAC.bp
  > rotate.pl -l 0 -L 4000 SEM/.BXE.semd.sac.bp



Note:
-----

the processing scripts make use of additional software
(assumed to be installed in a default directory /opt/seismo-util/):

- SAC
  http://www.iris.washington.edu/software/sac/

  dowload software package (IRIS DMC Software - Processing Programs)

  the scripts 'process_**.pl' make use of the provided binary `sac` and `saclst`.

- iaspei-tau
  http://www.iris.washington.edu/pub/programs/iaspei-tau/

  download software package (IRIS DMC Software - Processing Programs)
  and install it by running:
    > cd iaspei-tau/
    > ./build-package

  the scripts 'process_**.pl' call 'phtimes.csh' (for picking phases)
  which makes use of the provided binary `ttimes`.
  please modify the paths in the script 'phtimes.csh' accordingly.

- asc2sac
  source code provided in this directory or on https://geodynamics.org/portals/seismo/ inside
  the post-process utilities:
  https://geodynamics.org/portals/seismo/samples/seis_process-1.0.1.tar.gz

- convolve_stf
  source code provided in utils/lib or on https://geodynamics.org/portals/seismo/ inside
  the post-process utilities:
  https://geodynamics.org/portals/seismo/samples/seis_process-1.0.1.tar.gz



