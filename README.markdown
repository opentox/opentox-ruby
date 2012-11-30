opentox-ruby
============

Ruby wrapper for the [OpenTox](http://www.opentox.org) REST API 

Installation
------------

opentox-ruby depends on many third party programs and libraries, which makes the setup complicated and error prone. For this reason we recommend to use the installer from [opentox-install](http://github.com/opentox/install/tree/oldarch). If you want to install manually you can find the necessary steps in the installation scripts.

Quickstart
----------

This example shows how to create a lazar model and predict a compound, it assumes that you have access to a working installation of OpenTox services with corresponding settings in $HOME/.opentox/config. Run the following code in irb or from a ruby script:

    require 'rubygems'
    require 'opentox-ruby'

    # Authenticate
    subjectid = OpenTox::Authorization.authenticate(USER,PASSWORD) 

    # Upload a dataset
    training_dataset = OpenTox::Dataset.create_from_csv_file(TRAINING_DATASET, subjectid)

    # Create a prediction model
    model_uri = OpenTox::Algorithm::Lazar.new.run({:dataset_uri => training_dataset.uri, :subjectid => subjectid}).to_s
    lazar = OpenTox::Model::Lazar.find model_uri, subjectid
    
    # Predict a compound
    compound = OpenTox::Compound.from_smiles("c1ccccc1NN")
    prediction_uri = lazar.run(:compound_uri => compound.uri, :subjectid => subjectid)
    prediction = OpenTox::LazarPrediction.find(prediction_uri, subjectid)
    puts prediction.to_yaml

[API documentation](http://rdoc.info/gems/opentox-ruby/1.0.0/frames)
-------------------------------------------------------------------

Copyright
---------

Copyright (c) 2009-2012 Christoph Helma, Martin Guetlein, Micha Rautenberg, Andreas Maunz, David Vorgrimmler, Denis Gebele. See LICENSE for details.
