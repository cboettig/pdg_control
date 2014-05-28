---
layout: page
---


Abstract
--------

Optimal control entered mathematical ecology and evolution through applications in resource management and behavioral and evolutionary ecology. While traditional optimal control provides insights into problem structure and a set of analytic tools, it is most useful when combined with additional techniques that explicitly account for ill-resolved model structure, complex stochastics, and unusual or inflexible controls, which are commonly encountered in ecological systems. We will focus on three practical concerns that complicate straightforward applications of optimal control: model uncertainty, objective uncertainty, and control constraints. For model uncertainty, we consider scenarios where there are several competing models. Rather than applying optimal control to a model of dubious authority, a more reasonable strategy will reflect the uncertainty about the different models and update this through time as new information becomes available. For objective uncertainty, we consider cases where there is uncertainty about the correct objective function, such as managing a system for biological diversity with competing metrics of diversity. For control constraints, we consider fisheries management applications where optimal control solutions often represent functions that are either biologically or politically unfeasible, and we consider how to construct a principled approximation from within a class of pragmatically feasible controls. Using a problem-based approach, we will consider examples in fisheries management, marine spatial planning, life-history evolution, and ecosystem dynamics. Our team participants bring expertise from mathematics, ecology, control systems engineering, network theory, and statistics and are drawn from institutions in California, Pennsylvania, Hawaii, Tennessee, Florida, and Australia.


Current scripts
---------------

({{ site.pages | size }} entries)

{% for page in site.pages %}
- [{{ page.path }}](/{{site.repo}}{{ page.url }})
{% endfor %}

