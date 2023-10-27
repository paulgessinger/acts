%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THIS FILE IS AUTOGENERATED, DO NOT MANUALLY MODIFY IT!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Whitepapers

This is a collection of whitepapers documenting methods and approaches used
in ACTS.

{% for whp in config.whitepapers %}

## {{ whp.metadata.title }} ([GitHub]({{ whp.repository }}))

:::{admonition} Authors
:class: note
{% for author in whp.metadata.authors %}
- {{ author }}
{% endfor %}
:::


{{ whp.metadata.description }}

{% endfor %}