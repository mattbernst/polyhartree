The inchis.html provides a convenient one-page index of all species in the cccbdb collection.

The page generation logic is very flaky, hence the data is bundled here instead of fetched dynamically each time.

To retrieve an updated inchis.html, do:

curl http://cccbdb.nist.gov/inchi.asp > inchis.html

and repeat the attempt until the error message "Operation is not allowed" no longer shows up in the output.
