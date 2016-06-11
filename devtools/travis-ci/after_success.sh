#!/bin/bash
# Must be invoked with $PACKAGENAME

echo $TRAVIS_PULL_REQUEST $TRAVIS_BRANCH

if [ "$TRAVIS_PULL_REQUEST" = true ]; then
    echo "This is a pull request. No deployment will be done."; exit 0
fi


if [ "$TRAVIS_BRANCH" != "master" ]; then
    echo "No deployment on BRANCH='$TRAVIS_BRANCH'"; exit 0
fi


# Deploy to anaconda
conda install --yes anaconda-client jinja2 
pushd .
cd $HOME/miniconda/conda-bld
FILES=*/${PACKAGENAME}-dev-*.tar.bz2
for filename in $FILES; do
    #We're aiming to get this into omnia eventually but for now, dmobley channel
    anaconda -t $ANACONDA_TOKEN remove --force dmobley/${PACKAGENAME}-dev/${filename}
    anaconda -t $ANACONDA_TOKEN upload --force -u dmobley -p ${PACKAGENAME}-dev ${filename}
    #anaconda -t $ANACONDA_TOKEN remove --force ${ORGNAME}/${PACKAGENAME}-dev/${filename}
    #anaconda -t $ANACONDA_TOKEN upload --force -u ${ORGNAME} -p ${PACKAGENAME}-dev ${filename}
done
popd

