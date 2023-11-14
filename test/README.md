The test data exists in two archives in order to speed up testing:
../data/OrthoFinder_Test_Reference.tar.gz
../data/OrthoFinder_Test_Result_Reference.tar.gz

The first archive is for running with OrthoFinder and the second
archive is for running without OrthoFinder (i.e. it already contains
results from OrthoFinder). You set which file it is in the deploy.cfg,
in the test_data field, and you copy it to the ../test_local/refdata
directory.

It has to be unpacked. It should be unpacked, for testing purposes, in
../scripts/entrypoint.sh but I haven't figure out how to read
deploy.cfg in there yet. Also, the testing within the kb_orthofinder
app currently assumes that its the second archive, so you need to
force it to run, I might also do that in deploy.cfg.