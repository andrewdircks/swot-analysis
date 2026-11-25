- [x] what is SSH? SSHA?
    - [x] WGS84
    - [x] geoid
    - [x] "Geoid relief is about ±100 m: the ocean surface really does sit tens of metres higher over mass excesses"
    - [] geostrophic flow equations
    - [x] is SSHA an approximation of dynamic topography, that only exists becasue the geoid isn't granular enough? if SSHA was "perfect" would it equal exact dynamic topography (DT)?
        - NO ;; SSHA = DT - mean(DT)
        - absolute dynamic topography (ADT) = SSHA + mean(DT) -> THIS is the approximateion of dynamic topography; mean(DT) was independently estimate by special satellite missions, basically an estimate of longstanding currents
    - [x] if we want to use SSHA to measure persistent, large-scale currents (e.g. gulf stream), wouldn't mean sea surface (MSS) just converge to the steady state current, rather than isolate non-geoid caused elevation (of which the gulf stream is)?
        - use ADT instead
- [x] what is "SWOT pass"?
    - half revolution around the globe
- what is a cold eddy?
- [x] `SWOT_L2_LR_SSH_D`: what is this? what is the format of the data? where does it sit relative to other SWOT products? i.e. what are the levels?
- [x] what is `ssha_karin_2`?
    - direct SSH measurement from altimiter, adjusted for atmospheric + wave height
- [x] what is "crossover calibration"(`height_cor_xover`)?
    - don't really worry ... internal measurement adjustmenets

`swell_finder.py`
- [x] what is the "WindWave" granule? 
    - [x] what is the format of this data?
        - [x] what is the `swh_karin_qual` field?

`swell.py`
- [] what are the variables meanings? how were their values found?
    - [x] PASS ;; from "swell_finder.py + tile scan" ?
        - [x] "Pass 151 showed 30 m, but that's a coastal artifact"
            - [x] why does this coastal artifact show up
                - literally measuring mountains; still pass quality flag ...
            - [x] why are you now dropping 50km to coast?
                - arbitrary
        - [x] "a separate script that cut the pass into 40 km tiles and plotted each one's spectrum"
            - [x] why 40km?
                - want it large to capture many wavelengths; want it to fit in the swath
            - [x] what is the "spectrum"?
                - output of FT: wavelengths, orientations, amplitude
            - [x] why do we cut off the spectrum at 500m?
                - Nyquist limit: need >= 2 samples per wavelength (signal processing)
                - our sample rate is 250m
    - [] LAT_RANGE
        - [] "...about 35.5°N to 39.5°N, and only noise appeared north and south of that band"
            - [] noise? how could you tell?
        - [] "At 47.7°N, where SWH was largest, the waves were inside the storm and only noise showed."
            - [] "inside the storm": how do you know?
            - [] why would the largest SWH not be the clearest signal?
    - [] TILE_LAT 
    - [] TILE_KM
    - [] ZOOM_KM
    - [] BAND
- [] what is the Unsmoothed dataset?
    - [] what is the format of the data?
    - [] explain the 250m sample rate?

- [] FT refresh/deepdive
    - [x] vid https://www.youtube.com/watch?v=spUNpyF58BY
    - [] 2d FT
    - [] wave numbers (n/D)

