# Hints for Tutors: Tutorial 12 - Reaction-Ensemble, Constant-pH Method

## Learning objectives (physics):

After the tutorial, students should be able to 

* explain
    * difference between activity and concentration
    * the origin of the equilibrium constant K and the pK value
    * the origin of the pH value
    * basics of the reaction-esemble method and the constant pH method
       (implicit water and H+, insertion/deletion attempts, acceptance by probability based on energy difference)
* reason about the chemical interpretation of the neutralising B+ ion


## Learning objectives (ESPResSo)

After the tutorial, students should be able to 

* instance a reaction-ensemble constant-pH object und understand its parameters (temperature, exclusion radius)
* add reactions to the reaction-ensemble instance
* properly integrate the reactions into the integration/sampling scheme


## Points to mention throughout the tutorial

Make sure attendees understand 

* the (unphysical) concept of the exclusion radius
* how the exchange of particles works when performing a reaction
* that the H+ reservoir (and its particles) is not simulated explicitly but captured by the constant-pH method
* that the B+ ion is not neccessarily a H+ ion, but can depend on the system and pH-value
* that you can't blindly use the constant-pH-method for any pH, but it strongly depends on the exact system