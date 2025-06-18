## Database of Asteroid Models from Inversion Techniques (DAMIT)

DAMIT contains asteroid models that were derived using the light-curve inversion method developed by Kaasalainen & Torppa (2001) and Kaasalainen et al. (2001), combined with other inversion techniques in some cases. Each model in DAMIT references the original paper(s) where it was published. Please note that the models presented in DAMIT may differ from those published in the original papers. The main reason for this is a limited dataset used in the original publication and/or a narrow range of periods scanned during the inversion.

**Latest release preview available <a href="https://scraptechguy.github.io/DAMIT/index.html">here</a>.**

## My contribution to the project

#### Implemented light‐curve computation
+ Added `getLightCurve()` and other support functions to sample phase angles (γ) and compute corresponding fluxes from the 3D model, applying incidence/emergence cosines, shadowing masks, and Lambertian scattering.

+ Normalized the resulting γ and flux arrays for consistent plotting.

#### Integrated D3‐powered plotting
+ Built `drawLightCurve()` to render the phase–flux curve as an SVG line chart, complete with dynamic axes, padding, and responsive sizing alongside the existing asteroid projections.

+ Enhanced the UI with controls to scale the chart and to download light curve data as a two‐column TXT file or an SVG file.

<a href="https://astro.troja.mff.cuni.cz/projects/damit/">
<img src="https://github.com/user-attachments/assets/409c7306-51d5-472a-b39a-e8345e6e9fb9">
</a>

> © 2008–2025 <a href="https://astro.troja.mff.cuni.cz/">Astronomical Institute of the Charles University</a>, Josef Ďurech, Vojtěch Sidorin | Except where otherwise stated, content on this site is licensed under a <a href="https://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>. | Main contact: Josef Ďurech (durech@sirrah.troja.mff.cuni.cz)
