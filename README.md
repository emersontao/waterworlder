# **waterworlders code**

Starting fresh as of October 2026

**To-Do's:**
- Fix units issues in code IMMEDIATELY
- Use all previous _complete solutions to start; upload to VK

_Emerson_:
- comment code of biggrid.py
- add biggrid, windae_one, init weights to 
   `git add --force filepath` to add
- edit init.make_fname to allow for custom folder names
- Rework xpml, xp_evolution, etc. files


_Artem_:
- add automations to prevent biggrid.py from breaking while running on VK
- add an exception flag (separate from ramp_to, run_wind, etc)
- make it possible to save ramp_to and run_wind even if biggrid.py fails (separate file?)

- 3 total new/updated files:
  - biggrid, will take gravity as an input axis
  - biggrid_load_only, opens all files and constructs grid ONLY
  - the python script that creates a bash script to run n jobs in parallel for n different surface gravities

- Using Wind-AE 1.0: 
  - Create NEW span for the grids
    - cut off M_p at 25 M_e