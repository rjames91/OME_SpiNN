#top level SpiNNakEar vertex

#=============Needs to be able to overide partitioning and placement=============
# Needs to be able to produce the graph defined in model_launch_framework_binaural

#=============Needs to be initialised as a PyNN population=============
# pass in number of fibres as the n_neurons parameter and calculate all other parameters from there

#=============Needs to be able to send projections to other PyNN populations=============
# This means we have to send from the IHCAN vertices only
# on each core there are 2 IHCAN "neurons" so there are some options on how to proceed
# 1. Have a wrapper within the SpiNNakEar object that will create edges directly between each IHC machine vertex and
# connected neuron vertices when a from list connector is used. If we do this then we need to specify the number of atoms
# per vertex so the receiving neurons know how to interpret the incoming spikes - the downside of this approach is that
# the master population table size for any receiving neurons will be very large
# 2. An intermediate SpiNNakEar vertex is create that pools together the outputs from multiple IHCANs into one core
# This core can then transmit the received spikes as they arrive with its unique MV key | a converted neuron ID depending
# on which IHCAN the original spike came from.
# 3. We create a similar approach to 2. but incorporate the unique MV key + neuron IDs into the IHCAN c code. This
# will require some extra parameters to be passed down and I'm not sure how key generation will work but it seems
# like a more sensible solution than 2.

#=============Needs to be able to receive projections from other PyNN populations=============
# When a SpiNNakEar vertex receives a spike it must be routed to the corresponding DRNL instance
# If we are able to implement edge filtering and conn lut on DNRLs then I believe this could be straight forward
# The conn_lut in DTCM should be ok as this projection will always be from a single population consiting of <1000 neurons


#=============Needs to be recordable=============
# not essential to begin with but would like this feature in the future
#spinnakear is currently recordable on its own but this doesn't follow the standard PyNN recording methods.


#=============Questions================
# It seems that the routing key generation part of this is what I am most unsure about
# Essentially for PyNN integration the application vertex (SpiNNakEar) acts just like a wrapper
# hiding the interactions with individual vertices from the tools.
# e.g. when we have a projection to SpiNNakEar edges must be created between source neurons and the DRNL vertices
# there are fewer DRNL vertices than there are AN fibres (x10) so this must be made clear to a user
# To begin with I think forcing From list ocnnectors is the easiest implementation but further developements will need
# to enforce these differences when it comes to auto connection building.

#=============Unique parameters===========
#There will be some parameters to SpiNNakEar that will be unique, for example:
# number of ears
# audio stimulus/ sampling frequency
# each band frequency - this could be generated based on the number of fibres and the assumption of 10 fibres per band and a log dist.

