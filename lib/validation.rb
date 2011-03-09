module OpenTox
  class Validation
    include OpenTox
    
    # find validation, raises error if not found
    # @param [String] uri
    # @param [String,optional] subjectid
    # @return [OpenTox::Validation]
    def self.find( uri, subjectid=nil )
      val = Validation.new(uri)
      val.load_metadata( subjectid )
      val
    end
    
    # creates a validation object from crossvaldiation statistics, raise error if not found
    # (as crossvaldiation statistics are returned as an average valdidation over all folds)
    # @param [String] crossvalidation uri
    # @param [String,optional] subjectid
    # @return [OpenTox::Validation]
    def self.from_cv_statistics( crossvalidation_uri, subjectid=nil )
      find( File.join(crossvalidation_uri, 'statistics'),subjectid )
    end
    
    # loads metadata via yaml from validation object
    # fields (like for example the validated model) can be acces via validation.metadata[OT.model]
    def load_metadata( subjectid=nil )
      @metadata = YAML.load(OpenTox::RestClientWrapper.get(uri,{:subjectid => subjectid, :accept => "application/x-yaml"}))
    end
    
    # PENDING: creates summary as used for ToxCreate
    def summary
      if @metadata[OT.classificationStatistics]
        res = {
          :nr_predictions => @metadata[OT.numInstances] - @metadata[OT.numUnpredicted],
          :correct_predictions => @metadata[OT.classificationStatistics][OT.percentCorrect],
          :weighted_area_under_roc => @metadata[OT.classificationStatistics][OT.weightedAreaUnderRoc],
        }
        @metadata[OT.classificationStatistics][OT.classValueStatistics].each do |s|
          if s[OT.classValue].to_s=="true"
            res[:true_positives] = s[OT.numTruePositives]
            res[:false_positives] = s[OT.numFalsePositives]
            res[:true_negatives] = s[OT.numTrueNegatives]
            res[:false_negatives] = s[OT.numFalseNegatives]
            res[:sensitivity] = s[OT.truePositiveRate]
            res[:specificity] = s[OT.falsePositiveRate]
            break
          end
        end
        res
      elsif @metadata[OT.regressionStatistics]
        {
          :nr_predictions => @metadata[OT.numInstances] - @metadata[OT.numUnpredicted],
          :r_square => @metadata[OT.regressionStatistics][OT.rSquare],
          :root_mean_squared_error => @metadata[OT.regressionStatistics][OT.rootMeanSquaredError],
          :mean_absolute_error => @metadata[OT.regressionStatistics][OT.meanAbsoluteError],
        }
      end
    end
  end
  
  class Crossvalidation
    include OpenTox

    attr_reader :report
    
    # find crossvalidation, raises error if not found
    # @param [String] uri
    # @param [String,optional] subjectid
    # @return [OpenTox::Crossvalidation]
    def self.find( uri, subjectid=nil )
      cv = Crossvalidation.new(uri)
      cv.load_metadata( subjectid )
      cv
    end
		
    # creates a crossvalidations, waits until it finishes, may take some time
    # @param [Hash] params (required:algorithm_uri,dataset_uri,prediction_feature, optional:algorithm_params,num_folds(10),random_seed(1),stratified(false))
    # @param [String,optional] subjectid
    # @param [OpenTox::Task,optional] waiting_task (can be a OpenTox::Subtask as well), progress is updated accordingly
    # @return [OpenTox::Crossvalidation]
    def self.create( params, subjectid=nil, waiting_task=nil )
      params[:subjectid] = subjectid if subjectid
      uri = OpenTox::RestClientWrapper.post( File.join(CONFIG[:services]["opentox-validation"],"crossvalidation"),
        params,{:content_type => "text/uri-list"},waiting_task )
      Crossvalidation.new(uri)
    end

    # looks for report for this crossvalidation, creates a report if no report is found
    # @param [String,optional] subjectid
    # @param [OpenTox::Task,optional] waiting_task (can be a OpenTox::Subtask as well), progress is updated accordingly
    # @return [String] report uri
    def find_or_create_report( subjectid=nil, waiting_task=nil )
      @report = CrossvalidationReport.find_for_crossvalidation(@uri, subjectid) unless @report
      @report = CrossvalidationReport.create(@uri, subjectid, waiting_task) unless @report
      @report.uri
    end
    
    # loads metadata via yaml from crossvalidation object
    # fields (like for example the validations) can be acces via validation.metadata[OT.validation]
    def load_metadata( subjectid=nil )
      @metadata = YAML.load(OpenTox::RestClientWrapper.get(uri,{:subjectid => subjectid, :accept => "application/x-yaml"}))
    end
    
    # PENDING: creates summary as used for ToxCreate
    def summary( subjectid=nil )
      Validation.from_cv_statistics( @uri, subjectid ).summary
    end
  end
  
  class ValidationReport
    include OpenTox
    
    # finds ValidationReport for a particular validation
    # @param [String] crossvalidation uri 
    # @param [String,optional] subjectid
    # @return [OpenTox::ValidationReport] nil if no report found
    def self.find_for_validation( validation_uri, subjectid=nil )
      uris = RestClientWrapper.get(File.join(CONFIG[:services]["opentox-validation"],
        "/report/validation?validation="+validation_uri), {:subjectid => subjectid}).chomp.split("\n")
      uris.size==0 ? nil : ValidationReport.new(uris[-1])
    end
    
  end

  class CrossvalidationReport
    include OpenTox
    
    # finds CrossvalidationReport via uri, raises error if not found
    # @param [String] uri
    # @param [String,optional] subjectid
    # @return [OpenTox::CrossvalidationReport]
    def self.find( uri, subjectid=nil )
      # PENDING load report data?
      OpenTox::RestClientWrapper.get(uri,{:subjectid => subjectid})
      CrossvalidationReport.new(uri)
    end
    
    # finds CrossvalidationReport for a particular crossvalidation
    # @param [String] crossvalidation uri 
    # @param [String,optional] subjectid
    # @return [OpenTox::CrossvalidationReport] nil if no report found
    def self.find_for_crossvalidation( crossvalidation_uri, subjectid=nil )
      uris = RestClientWrapper.get(File.join(CONFIG[:services]["opentox-validation"],
        "/report/crossvalidation?crossvalidation="+crossvalidation_uri), {:subjectid => subjectid}).chomp.split("\n")
      uris.size==0 ? nil : CrossvalidationReport.new(uris[-1])
    end
    
    # creates a crossvalidation report via crossvalidation
    # @param [String] crossvalidation uri 
    # @param [String,optional] subjectid
    # @param [OpenTox::Task,optional] waiting_task (can be a OpenTox::Subtask as well), progress is updated accordingly
    # @return [OpenTox::CrossvalidationReport]
    def self.create( crossvalidation_uri, subjectid=nil, waiting_task=nil )
      uri = RestClientWrapper.post(File.join(CONFIG[:services]["opentox-validation"],"/report/crossvalidation"),
        { :validation_uris => crossvalidation_uri, :subjectid => subjectid }, {}, waiting_task )
      CrossvalidationReport.new(uri)
    end
  end
  
  class QMRFReport
    include OpenTox
    
    # finds QMRFReport, raises Error if not found
    # @param [String] uri
    # @param [String,optional] subjectid
    # @return [OpenTox::QMRFReport]
    def self.find( uri, subjectid=nil )
      # PENDING load crossvalidation data?
      OpenTox::RestClientWrapper.get(uri,{:subjectid => subjectid})
      QMRFReport.new(uri)
    end
    
    # finds QMRF report for a particular model
    # @param [String] model_uri 
    # @param [String,optional] subjectid
    # @return [OpenTox::QMRFReport] nil if no report found
    def self.find_for_model( model_uri, subjectid=nil )
      uris = RestClientWrapper.get(File.join(CONFIG[:services]["opentox-validation"],
        "/reach_report/qmrf?model="+model_uri), {:subjectid => subjectid}).chomp.split("\n")
      uris.size==0 ? nil : QMRFReport.new(uris[-1])
    end
    
    # creates a qmrf report via model
    # @param [String] model_uri 
    # @param [String,optional] subjectid
    # @param [OpenTox::Task,optional] waiting_task (can be a OpenTox::Subtask as well), progress is updated accordingly
    # @return [OpenTox::QMRFReport]
    def self.create( model_uri, subjectid=nil, waiting_task=nil )
      uri = RestClientWrapper.post(File.join(CONFIG[:services]["opentox-validation"],"/reach_report/qmrf"), 
        { :model_uri => model_uri, :subjectid => subjectid }, {}, waiting_task )
      QMRFReport.new(uri)
    end
  end
  
end

