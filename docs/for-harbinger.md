# For the Harbinger team

The part of this repo that matters for a battery-electric work truck is the coupling of pack temperature to sunlight, duty cycle and derating:

- A two-node transient thermal model driven by real hourly weather (Dallas, 2025), with solar absorptance and IR emissivity of the skin as explicit inputs.
- Arrhenius calendar-aging and BMS derating, integrated over a full year.
- A film-on-substrate optics model that says which coatings or laminates actually lower solar gain. ([design map](figures/laminate_design_map.png))

**A 3-week project I could do for you, remote:** swap the grid-dispatch duty cycle for a delivery-route duty cycle (parked in sun, then driving), and estimate how body color/coating and parking strategy change pack temperature, derating and aging across hot-climate depots.

— Aneesh Karavadi · aneesh.karavadi@gmail.com
