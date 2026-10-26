# NGC 3201 orbit inputs

The [source record](ngc3201-source.json) preserves one row from the Baumgardt database's **/newdata/orbits.html** table, captured 25 October 2026. The [unaltered row extract](ngc3201-row.html) retains the native values and provenance attributes; the complete original HTML is cached locally under `raw/` and identified by SHA-256. It is not mixed with the root-path table, which contains different values.

The record includes observed sky coordinates, distance, radial velocity and mean proper motion, plus model-derived Galactocentric Cartesian coordinates and velocities. The latter retain their published frame and rounding. No complete six-dimensional covariance or complete solar-frame prescription is supplied by this table. Proper-motion correlation and quoted errors are retained; they are not replaced by invented independent errors on the six Cartesian components.

Changing distance in the defined sensitivity cases adds the corresponding heliocentric displacement and transverse-velocity difference to the catalog Cartesian anchor. It does not refit the solar frame. The selected host's distance is a separate conditional case, not a joint posterior formed by multiplying overlapping datasets.

The Galactic potential uses the numerical parameters of **Model I, Table 1**, in [Irrgang et al. (2013)](https://doi.org/10.1051/0004-6361/201220540), with the exact potential unit 100 (km/s)^2. This avoids converting through the rounded mass unit printed in the article. The catalog itself names Irrgang et al. for its orbit calculations; exact equality of integration, averaging and source versions is not assumed.

Original data and publication rights remain with their sources. The repository license applies to the original extraction/analysis code and documentation.
