module OpenTox
  
  class Dataset 

    attr_accessor :uri, :title, :creator, :data, :features, :compounds

    def initialize( owl=nil )
      @data = {}
      @features = []
      @compounds = []
      
      # creates dataset object from Opentox::Owl object
      # use Dataset.find( <uri> ) to load dataset from rdf-supporting datasetservice
      # note: does not load all feature values, as this is time consuming
      if owl
        raise "invalid param" unless owl.is_a?(OpenTox::Owl)
        @title = owl.get("title")
        @creator = owl.get("creator")
        @uri = owl.uri
        # when loading a dataset from owl, only compound- and feature-uris are loaded 
        owl.load_dataset(@compounds, @features)
        # all features are marked as dirty
        # as soon as a feature-value is requested all values for this feature are loaded from the rdf
        @dirty_features = @features.dclone
        @owl = owl
      end
    end
  
    # Find and a dataset and load all data from the dataset service
    # @param [String] uri Dataset URI
    # @return [OpenTox::Dataset] Dataset object
    def self.find(uri, accept_header=nil) 
    
      unless accept_header
        if (@@config[:yaml_hosts].include?(URI.parse(uri).host))
          accept_header = 'application/x-yaml'
        else
          accept_header = "application/rdf+xml"
        end
      end
      
      case accept_header
      when "application/x-yaml"
        d = YAML.load RestClientWrapper.get(uri.to_s.strip, :accept => 'application/x-yaml').to_s 
        d.uri = uri unless d.uri
      when "application/rdf+xml"
        owl = OpenTox::Owl.from_uri(uri.to_s.strip, "Dataset")
        d = Dataset.new(owl)
      else
        raise "cannot get datset with accept header: "+accept_header.to_s
      end
      d
    end
    
    # converts a dataset represented in owl to yaml
    # (uses a temporary dataset)
    # note: to_yaml is overwritten, loads complete owl dataset values 
    def self.owl_to_yaml( owl_data, uri)
      owl = OpenTox::Owl.from_data(owl_data, uri, "Dataset")
      d = Dataset.new(owl)
      d.to_yaml
    end
    
    # creates a new dataset, using only those compounsd specified in new_compounds
    # returns uri of new dataset
    def create_new_dataset( new_compounds, new_features, new_title, new_creator )
      
      LOGGER.debug "create new dataset with "+new_compounds.size.to_s+"/"+compounds.size.to_s+" compounds"
      raise "no new compounds selected" unless new_compounds and new_compounds.size>0
      
      # load require features 
      if ((defined? @dirty_features) && (@dirty_features & new_features).size > 0)
        (@dirty_features & new_features).each{|f| load_feature_values(f)}
      end
      
      dataset = OpenTox::Dataset.new
      dataset.title = new_title
      dataset.creator = new_creator
      dataset.features = new_features
      dataset.compounds = new_compounds
      
      # Copy dataset data for compounds and features
      # PENDING: why storing feature values in an array? 
      new_compounds.each do |c|
        data_c = []
        raise "no data for compound '"+c.to_s+"'" if @data[c]==nil
        @data[c].each do |d|
          m = {}
          new_features.each do |f|
            m[f] = d[f]
          end
          data_c << m 
        end
        dataset.data[c] = data_c
      end
      return dataset.save
    end
    
    # returns classification value
    def get_predicted_class(compound, feature)
      v = get_value(compound, feature)
      if v.is_a?(Hash)
        k = v.keys.grep(/classification/).first
        unless k.empty?
        #if v.has_key?(:classification)
          return v[k]
        else
          return "no classification key"
        end
      elsif v.is_a?(Array)
        raise "predicted class value is an array\n"+
          "value "+v.to_s+"\n"+
          "value-class "+v.class.to_s+"\n"+
          "dataset "+@uri.to_s+"\n"+
          "compound "+compound.to_s+"\n"+
          "feature "+feature.to_s+"\n"
      else
        return v
      end
    end
    
    # returns regression value
    def get_predicted_regression(compound, feature)
      v = get_value(compound, feature)
      if v.is_a?(Hash)
        k = v.keys.grep(/regression/).first
        unless k.empty?
          return v[k]
        else
          return "no regression key"
        end
      elsif v.is_a?(Array)
        raise "predicted regression value is an array\n"+
          "value "+v.to_s+"\n"+
          "value-class "+v.class.to_s+"\n"+
          "dataset "+@uri.to_s+"\n"+
          "compound "+compound.to_s+"\n"+
          "feature "+feature.to_s+"\n"
      else
        return v
      end
    end
    
    # returns prediction confidence if available
    def get_prediction_confidence(compound, feature)
      v = get_value(compound, feature)
      if v.is_a?(Hash)
        k = v.keys.grep(/confidence/).first
        unless k.empty?
        #if v.has_key?(:confidence)
          return v[k].abs
          #return v["http://ot-dev.in-silico.ch/model/lazar#confidence"].abs
        else
          # PENDING: return nil isntead of raising an exception
          raise "no confidence key"
        end
      else
        LOGGER.warn "no confidence for compound: "+compound.to_s+", feature: "+feature.to_s
        return 1
      end
    end
    
    # return compound-feature value
    def get_value(compound, feature)
      if (defined? @dirty_features) && @dirty_features.include?(feature)
        load_feature_values(feature)
      end
      
      v = @data[compound]
      return nil if v == nil # missing values for all features
      if v.is_a?(Array)
        # PENDING: why using an array here?
        v.each do |e|
          if e.is_a?(Hash)
            if e.has_key?(feature)
              return e[feature]
            end
          else
            raise "invalid internal value type"
          end
        end
        return nil #missing value
      else
        raise "value is not an array\n"+
              "value "+v.to_s+"\n"+
              "value-class "+v.class.to_s+"\n"+
              "dataset "+@uri.to_s+"\n"+
              "compound "+compound.to_s+"\n"+
              "feature "+feature.to_s+"\n"
      end
    end

    # loads specified feature and removes dirty-flag, loads all features if feature is nil
    def load_feature_values(feature=nil)
      if feature
        raise "feature already loaded" unless @dirty_features.include?(feature)
        @owl.load_dataset_feature_values(@compounds, @data, [feature])
        @dirty_features.delete(feature)
      else
        @data = {} unless @data
        @owl.load_dataset_feature_values(@compounds, @data, @dirty_features)
        @dirty_features.clear
      end
    end
    
    # overwrite to yaml:
    # in case dataset is loaded from owl:
    # * load all values 
    def to_yaml
      # loads all features  
      if ((defined? @dirty_features) && @dirty_features.size > 0)
        load_feature_values
      end
      super
    end
    
    # * remove @owl from yaml, not necessary
    def to_yaml_properties
      super - ["@owl"]
    end

    # saves (changes) as new dataset in dataset service
    # returns uri
    # uses to yaml method (which is overwritten)
    def save
      OpenTox::RestClientWrapper.post(@@config[:services]["opentox-dataset"],{:content_type =>  "application/x-yaml"},self.to_yaml).strip   
    end
  end
end
