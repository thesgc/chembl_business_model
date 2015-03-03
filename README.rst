.. image:: https://badge.waffle.io/chembl/chembl_business_model.png?label=ready&title=Ready 
 :target: https://waffle.io/chembl/chembl_business_model
 :alt: 'Stories in Ready'
chembl_business_model
======

.. image:: https://pypip.in/v/chembl_business_model/badge.png
    :target: https://crate.io/packages/chembl_business_model/
    :alt: Latest PyPI version

.. image:: https://pypip.in/d/chembl_business_model/badge.png
    :target: https://crate.io/packages/chembl_business_model/
    :alt: Number of PyPI downloads

This is chembl_business_model package developed at Chembl group, EMBL-EBI, Cambridge, UK.

This package provides a django ORM model decorated with additional information about various aspects like searching and api exposing.
It also adds some application specific logic like what to do when the compound is saved.

The package uses 'chembl_core_model' as its base and declares new model by proxy inheritance.
This model is required when you want to create/update (write/modify) objects to database through ORM.
Core model intended only for read, with methods from this module you can also write.

This module defines also some application specific models (in meta.py file) which can be kept in separate db schema.
