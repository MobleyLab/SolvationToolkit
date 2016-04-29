This is a recipe for building the current development package into a conda
binary.

The installation on travis-ci is done by building the conda package, installing
it, running the tests, and then if successful (the idea would be to proceed by) pushing the package to binstar. The binstar auth token is an encrypted environment variable generated using:

binstar auth -n openmoltools-travis -o omnia --max-age 22896000 -c --scopes api:write

and then saved in the environment variable BINSTAR_TOKEN.

Note that most of the above is from openmoltools and we are currently working on adapting for solvationtoolkit.


