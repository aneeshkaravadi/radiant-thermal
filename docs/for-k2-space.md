# For the K2 Space team

High-power satellites live or die by radiator area. `radiator.py` sizes it from first principles:

$$A = \frac{Q}{\varepsilon\sigma T^4 - [\alpha S\cos\theta + F_E(\varepsilon q_{IR} + \alpha a S)]}$$

At 300 K in a LEO hot case, an optical solar reflector needs 4.1 m² per kW and white paint 5.3 m². Black paint cannot reject heat. ([figure](figures/radiator_area.png))

The same α/ε balance runs the rest of the repo, from rooftop radiative cooling to battery enclosures driven by real weather.

**A 3-week project I could do for you, remote:** an orbit-resolved radiator sizing tool (beta angle, eclipse, attitude), with end-of-life coating degradation and a mass estimate per kW, for early bus trades.

— Aneesh Karavadi · aneesh.karavadi@gmail.com
