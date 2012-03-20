@@cactus_uri="http://cactus.nci.nih.gov/chemical/structure/"
@@ambit_uri="http://ambit.uni-plovdiv.bg:8080/ambit2/depict/cdk?search="

module OpenTox

  require "rexml/document"
  # Ruby wrapper for OpenTox Compound Webservices (http://opentox.org/dev/apis/api-1.2/structure).
	class Compound 

    include OpenTox

		attr_accessor :inchi, :uri

		# Create compound with optional uri
    # @example
    #   compound = Compound.new("http://webservices.in-silico.ch/compound/InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"")
    # @param [optional, String] uri Compound URI
    # @return [Compound] Compound
		def initialize(uri=nil)
      @uri = uri
      case @uri
      when /InChI/ # shortcut for IST services
        @inchi = @uri.sub(/^.*InChI/, 'InChI')
      else
        @inchi = RestClientWrapper.get(@uri, :accept => 'chemical/x-inchi').to_s.chomp if @uri
      end
      
      if @uri and @inchi.to_s.size==0
        LOGGER.warn "REMOVE ABMIT HACK: no inchi for compound "+@uri.to_s+", load via smiles"
        @inchi = Compound.smiles2inchi(Compound.smiles(@uri))
      end
    end
    
    # request smiles from compound service via accept header
    # @return smiles as string
    def self.smiles(uri)
      RestClientWrapper.get(uri, :accept => 'chemical/x-daylight-smiles').to_s.chomp
    end

    # Create a compound from smiles string
    # @example
    #   compound = Compound.from_smiles("c1ccccc1")
    # @param [String] smiles Smiles string
    # @return [Compound] Compound
    def self.from_smiles(smiles)
      c = Compound.new
      c.inchi = Compound.smiles2inchi(smiles)
      c.uri = File.join(CONFIG[:services]["opentox-compound"],URI.escape(c.inchi))
      c
    end

    # Create a compound from inchi string
    # @param [String] smiles InChI string
    # @return [Compound] Compound
    def self.from_inchi(inchi)
      c = Compound.new
      c.inchi = inchi
      c.uri = File.join(CONFIG[:services]["opentox-compound"],URI.escape(c.inchi))
      c
    end

    # Create a compound from sdf string
    # @param [String] smiles SDF string
    # @return [Compound] Compound
    def self.from_sdf(sdf)
      c = Compound.new
      c.inchi = Compound.sdf2inchi(sdf)
      c.uri = File.join(CONFIG[:services]["opentox-compound"],URI.escape(c.inchi))
      c
    end

    # Create a compound from name. Relies on an external service for name lookups.
    # @example
    #   compound = Compound.from_name("Benzene")
    # @param [String] name name can be also an InChI/InChiKey, CAS number, etc
    # @return [Compound] Compound
    def self.from_name(name)
      c = Compound.new
      # paranoid URI encoding to keep SMILES charges and brackets
      c.inchi = RestClientWrapper.get("#{@@cactus_uri}#{URI.encode(name, Regexp.new("[^#{URI::PATTERN::UNRESERVED}]"))}/stdinchi").to_s.chomp
      c.uri = File.join(CONFIG[:services]["opentox-compound"],URI.escape(c.inchi))
      c
    end

		# Get InChI
    # @return [String] InChI string
		def to_inchi
      @inchi
		end

		# Get (canonical) smiles
    # @return [String] Smiles string
		def to_smiles
			Compound.obconversion(@inchi,'inchi','can')
		end

    # Get sdf
    # @return [String] SDF string
		def to_sdf
			Compound.obconversion(@inchi,'inchi','sdf')
		end

    # Get gif image
    # @return [image/gif] Image data
		def to_gif
			RestClientWrapper.get("#{@@cactus_uri}#{@inchi}/image")
		end

    # Get png image
    # @example
    #   image = compound.to_png
    # @return [image/png] Image data
		def to_png
      RestClientWrapper.get(File.join @uri, "image")
		end

    # Get URI of compound image
    # @return [String] Compound image URI
		def to_image_uri
      File.join @uri, "image"
		end

    # Get all known compound names. Relies on an external service for name lookups.
    # @example
    #   names = compound.to_names
    # @return [String] Compound names
		def to_names
      begin
        RestClientWrapper.get("#{@@cactus_uri}#{@inchi}/names").split("\n")
      rescue
        "not available"
      end
		end
    
    
    # Get all known compound names sorted by classification. Relies on an external service for name lookups.
    # @example
    #   names = compound.to_names_hash
    # @return [Hash] Classification => Name Array
		def to_names_hash
      begin
        xml = RestClientWrapper.get("#{@@cactus_uri}#{@inchi}/names/xml")
        xmldoc = REXML::Document.new(xml)
        data = {}
        
        xmldoc.root.elements[1].elements.each{|e|
          if data.has_key?(e.attribute("classification").value) == false
             data[e.attribute("classification").value] = [e.text]
          else
             data[e.attribute("classification").value].push(e.text)
          end
        }
        data
      rescue
        "not available"
      end
		end

    # Get all known compound names sorted by classification. Relies on an external service for name lookups.
    # @example
    #   names = compound.to_names_hash
    # @return [Hash] Classification => Name Array
    def to_ambit_names_hash
      begin
        ds = OpenTox::Dataset.new
        ds.save
        ds.load_rdfxml(RestClientWrapper.get("http://apps.ideaconsult.net:8080/ambit2/query/compound/search/names?type=smiles&property=&search=#{@inchi}"))
        ds.save
        ds.uri
      rescue
        "not available"
      end
    end


		# Match a smarts string
    # @example
    #   compound = Compound.from_name("Benzene")
    #   compound.match?("cN") # returns false
    # @param [String] smarts Smarts string
		def match?(smarts)
			obconversion = OpenBabel::OBConversion.new
			obmol = OpenBabel::OBMol.new
			obconversion.set_in_format('inchi')
			obconversion.read_string(obmol,@inchi) 
			smarts_pattern = OpenBabel::OBSmartsPattern.new
			smarts_pattern.init(smarts)
			smarts_pattern.match(obmol)
		end

		# Match an array of smarts strings, returns array with matching smarts
    # @example
    #   compound = Compound.from_name("Benzene")
    #   compound.match(['cc','cN']) # returns ['cc']
    # @param [Array] smarts_array Array with Smarts strings
    # @return [Array] Array with matching Smarts strings
		def match(smarts_array)
      # avoid recreation of OpenBabel objects
			obconversion = OpenBabel::OBConversion.new
			obmol = OpenBabel::OBMol.new
			obconversion.set_in_format('inchi')
			obconversion.read_string(obmol,@inchi) 
			smarts_pattern = OpenBabel::OBSmartsPattern.new
			smarts_array.collect do |smarts|
        smarts_pattern.init(smarts)
        smarts if smarts_pattern.match(obmol)
      end.compact
      #smarts_array.collect { |s| s if match?(s)}.compact
		end

    # Match_hits an array of smarts strings, returns hash with matching smarts as key and number of non-unique hits as value 
    # @example
    #   compound = Compound.from_name("Benzene")
    #   compound.match(['cc','cN']) # returns ['cc']
    # @param [Array] smarts_array Array with Smarts strings
    # @return [Hash] Hash with matching smarts as key and number of non-unique hits as value
		def match_hits(smarts_array)
      # avoid recreation of OpenBabel objects
			obconversion = OpenBabel::OBConversion.new
			obmol = OpenBabel::OBMol.new
			obconversion.set_in_format('inchi')
			obconversion.read_string(obmol,@inchi) 
			smarts_pattern = OpenBabel::OBSmartsPattern.new
			smarts_hits = {}
      #LOGGER.debug "dv ----------- obmol  #{Compound.new(@inchi).to_smiles}"
      smarts_array.collect do |smarts|
        #LOGGER.debug "dv ----------- all smarts  #{smarts}"
        smarts_pattern.init(smarts)
        if smarts_pattern.match(obmol)
          hits = smarts_pattern.get_map_list
          smarts_hits[smarts] = hits.size 
        end
      end
      #LOGGER.debug "dv ----------- smarts => hits #{smarts_hits}"
      return smarts_hits
      #smarts_array.collect { |s| s if match?(s)}.compact
		end
    
    # Lookup numerical values, returns hash with feature name as key and value as value 
    # @param [Array] Array of feature names
    # @param [String] Feature dataset uri
    # @param [String] Comma separated pc types
    # @return [Hash] Hash with feature name as key and value as value
		def lookup(feature_array,feature_dataset_uri,pc_type)
      ds = OpenTox::Dataset.find(feature_dataset_uri)

      #entry = ds.data_entries[self.uri]
      entry = nil 
      ds.data_entries.each { |c_uri, values| 
        if c_uri.split('/compound/').last == self.to_inchi
          entry = ds.data_entries[self.uri]
          break
        end
      }

      LOGGER.debug "#{entry.size} entries in feature ds for query." unless entry.nil?

      if entry.nil?
        temp_ds = OpenTox::Dataset.create; temp_ds.add_compound(self.uri)
        uri = RestClientWrapper.post(temp_ds.save + "/pcdesc", {:pc_type => pc_type})
        ds = OpenTox::Dataset.find(uri)
        entry = ds.data_entries[self.uri]
        ds.delete
        temp_ds.delete
      end

      features = entry.keys
      features.each { |feature| 
        new_feature = File.join(feature_dataset_uri, "feature", feature.split("/").last) 
        entry[new_feature] = entry[feature].flatten.first.to_f # see algorithm/lazar.rb:182, to_f because feature type detection doesn't work w 1 entry
        entry.delete(feature) unless feature == new_feature # e.g. when loading from ambit
      }
      #res = feature_array.collect {|v| entry[v]}
      entry
		end


    # Get URI of compound image with highlighted fragments
    #
    # @param [Array] activating Array with activating Smarts strings
    # @param [Array] deactivating Array with deactivating Smarts strings
    # @return [String] URI for compound image with highlighted fragments
    def matching_smarts_image_uri(activating, deactivating)
      activating_smarts = "\"#{activating.collect{|smarts| CGI.escape(smarts)}.join("\"/\"")}\""
      deactivating_smarts = "\"#{deactivating.collect{|smarts| CGI.escape(smarts)}.join("\"/\"")}\""
      File.join @uri, "smarts/activating", activating_smarts, "deactivating", deactivating_smarts
    end


    private

    # Convert sdf to inchi
		def self.sdf2inchi(sdf)
			Compound.obconversion(sdf,'sdf','inchi')
		end

    # Convert smiles to inchi
		def self.smiles2inchi(smiles)
			Compound.obconversion(smiles,'smi','inchi')
		end

    # Convert smiles to canonical smiles
		def self.smiles2cansmi(smiles)
			Compound.obconversion(smiles,'smi','can')
		end

    # Convert identifier from OpenBabel input_format to OpenBabel output_format
		def self.obconversion(identifier,input_format,output_format)
			obconversion = OpenBabel::OBConversion.new
			obmol = OpenBabel::OBMol.new
			obconversion.set_in_and_out_formats input_format, output_format
			obconversion.read_string obmol, identifier
			case output_format
			when /smi|can|inchi/
				obconversion.write_string(obmol).gsub(/\s/,'').chomp
			else
				obconversion.write_string(obmol)
			end
		end
	end
end
