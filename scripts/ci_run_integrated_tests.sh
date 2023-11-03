# Downloads the prebuit binary of GEOS and run the integrated tests.
# Then upload the results to a bucket such that those results can be further used.

# TODO we could unpack geos in a clean ubuntu image.

echo "hello there!"

echo ${SHORT_COMMIT}
# Overwrite the TPLs as well (disputable)
# curl -fSL https://storage.googleapis.com/geosx/ubuntu22.04-gcc11/GEOSX-and-TPL-${SHORT_COMMIT}.tar.gz | tar --directory=${GEOSX_TPL_DIR} --strip-components=1 --skip-old-files -xvzf -
curl -fSL https://storage.googleapis.com/geosx/ubuntu22.04-gcc11/GEOSX-and-TPL-${SHORT_COMMIT}.tar.gz | tar --directory=${GEOSX_TPL_DIR} --strip-components=1 -xvzf -
