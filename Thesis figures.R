#a systems diagram/ToC/feedback loop diagram for LNRS
install.packages("DiagrammeR")
library(DiagrammeR)

grViz("

digraph LNRS {

graph [layout = dot]

node [
shape = box
style = filled
fillcolor = black
fontcolor = white
fontsize = 12
]

Drivers -> LNRS
LNRS -> Habitat
Habitat -> Connectivity
Connectivity -> Biodiversity
Biodiversity -> Services
Services -> Wellbeing
Wellbeing -> Support
Support -> Investment
Investment -> LNRS

Drivers [label='Drivers of Decline']
LNRS [label='LNRS Actions']
Habitat [label='Habitat Creation & Restoration']
Connectivity [label='Ecological Connectivity']
Biodiversity [label='Biodiversity Recovery']
Services [label='Ecosystem Services']
Wellbeing [label='Human Wellbeing']
Support [label='Public & Political Support']
Investment [label='Funding & Partnerships']

}
")

#a better version thesis standard
library(DiagrammeR)

grViz("

digraph LNRS {

graph [
layout = dot
rankdir = TB
]

node [
shape = rectangle
style = 'rounded,filled'
fillcolor = honeydew
color = darkgreen
fontname = Helvetica
fontsize = 11
]

Drivers      [label='Environmental Pressures\nHabitat loss\nPollution\nClimate change']
Actions      [label='LNRS Actions']
Networks     [label='Habitat Restoration\nand Connectivity']
Biodiversity [label='Biodiversity Recovery']
Services     [label='Ecosystem Services']
Wellbeing    [label='Human Wellbeing']
Support      [label='Public & Political Support']
Investment   [label='Funding and Partnerships']

Drivers -> Actions
Actions -> Networks
Networks -> Biodiversity
Biodiversity -> Services
Services -> Wellbeing
Wellbeing -> Support
Support -> Investment
Investment -> Actions [label='reinforcing loop']

}
")

#three loops systems diagram allowing us to discuss LNRS 
#as a socioecological system rather than a simple intervention model

library(DiagrammeR)

grViz("

digraph LNRS {

graph [
  rankdir = LR
  overlap = false
]

node [
  shape = box
  style = 'rounded,filled'
  fillcolor = white
  fontname = Helvetica
  fontsize = 32
]

# =====================
# Nodes
# =====================

Actions [
 label='LNRS Actions\nRestoration, creation,\nmanagement'
]

Networks [
 label='Habitat Networks\nMore, Bigger, Better,\nMore Joined-Up'
]

Biodiversity [
 label='Biodiversity\nRecovery'
]

Services [
 label='Ecosystem\nServices'
]

Benefits [
 label='Human & Economic\nBenefits'
]

Support [
 label='Public & Political\nSupport'
]

Investment [
 label='Funding &\nPartnerships'
]

Monitoring [
 label='Monitoring'
]

Evidence [
 label='Evidence &\nLearning'
]

Adaptation [
 label='Adaptive\nManagement'
]

# =====================
# R1 Ecological Recovery
# =====================

Actions -> Networks [
 color='forestgreen'
 penwidth=4
]

Networks -> Biodiversity [
 color='forestgreen'
 penwidth=4
]

Biodiversity -> Services [
 color='forestgreen'
 penwidth=4
]

Services -> Networks [
 color='forestgreen'
 penwidth=4
 label='R1'
 fontcolor='forestgreen'
]

# =====================
# R2 Social Support
# =====================

Services -> Benefits [
 color='steelblue'
 penwidth=4
]

Benefits -> Support [
 color='steelblue'
 penwidth=4
]

Support -> Investment [
 color='steelblue'
 penwidth=4
]

Investment -> Actions [
 color='steelblue'
 penwidth=4
 label='R2'
 fontcolor='steelblue'
]

# =====================
# R3 Adaptive Governance
# =====================

Monitoring -> Evidence [
 color='darkorange'
 penwidth=4
]

Evidence -> Adaptation [
 color='darkorange'
 penwidth=4
]

Adaptation -> Actions [
 color='darkorange'
 penwidth=4
]

Actions -> Monitoring [
 color='darkorange'
 penwidth=4
 label='R3'
 fontcolor='darkorange'
]

# =====================
# Loop Labels
# =====================

Loop1 [
  shape=plaintext
  fontcolor='forestgreen'
  label='R1 Ecological Recovery'
]

Loop2 [
  shape=plaintext
  fontcolor='steelblue'
  label='R2 Social Support'
]

Loop3 [
  shape=plaintext
  fontcolor='darkorange'
  label='R3 Adaptive Governance'
]

}
")
#The Bedfordshire LNRS is underpinned by a strong ecological 
#theory of change (R1), contains an implicit social 
#legitimacy and investment feedback (R2), but lacks a 
#clearly articulated adaptive governance loop (R3), 
#limiting its operationalisation of resilience as a 
#measurable and continuously learned property of the system.