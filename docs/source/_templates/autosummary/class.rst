{{ fullname | escape | underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
   :members:
   :inherited-members:

   {% block methods %}
   {% set visible_methods = methods | reject("in", inherited_members) | reject("eq", "__init__") | list %}
   {% if visible_methods %}
   .. rubric:: Methods

   .. autosummary::
      :nosignatures:
   {% for item in visible_methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}
