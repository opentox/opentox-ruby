module OpenTox
	class Crossvalidation
    include OpenTox

		attr_reader :report
    
    # find crossvalidation, raises error if not found
    # @param [String] uri
    # @param [String,optional] subjectid
    # @return [OpenTox::Crossvalidation]
    def self.find( uri, subjectid=nil )
      # PENDING load crossvalidation data?
      OpenTox::RestClientWrapper.get(uri,{:subjectid => subjectid})
      Crossvalidation.new(uri)
    end
		
    # creates a crossvalidations, waits until it finishes, may take some time
    # @param [Hash] params
    # @param [String,optional] subjectid
    # @param [OpenTox::Task,optional] waiting_task (can be a OpenTox::Subtask as well), progress is updated accordingly
    # @return [OpenTox::Crossvalidation]
    def self.create( params, subjectid=nil, waiting_task=nil )
      params[:uri] = File.join(CONFIG[:services]['opentox-validation'], "crossvalidation")
      params[:num_folds] = 10 unless params[:num_folds]
=begin
<<<<<<< HEAD:lib/validation.rb
		 	params[:random_seed] = 2 unless params[:random_seed]
		 	params[:stratified] = false unless params[:stratified]
      uri = OpenTox::RestClientWrapper.post(File.join(CONFIG[:services]["opentox-validation"],"/crossvalidation"),params,nil,false)
			OpenTox::Validation.new(uri)
		end

    def create_report(subjectid=nil)
			@report_uri = OpenTox::RestClientWrapper.post(File.join(CONFIG[:services]["opentox-validation"],"/report/crossvalidation"), {:validation_uris => @uri, :subjectid => subjectid}).to_s
      @report_uri
    end

    def create_qmrf_report(subjectid=nil)
      @qmrf_report_uri = OpenTox::RestClientWrapper.post(File.join(CONFIG[:services]["opentox-validation"],"/reach_report/qmrf"), {:model_uri => @uri, :subjectid => subjectid}).to_s
      @qmrf_report_uri
    end

    def summary(type, subjectid=nil)
      v = YAML.load OpenTox::RestClientWrapper.get(File.join(@uri, 'statistics'),{:accept => "application/x-yaml", :subjectid => subjectid}).to_s

      case type
      when "classification"
        tp=0; tn=0; fp=0; fn=0; n=0
        v[:classification_statistics][:confusion_matrix][:confusion_matrix_cell].each do |cell|
          if cell[:confusion_matrix_predicted] == "true" and cell[:confusion_matrix_actual] == "true"
            tp = cell[:confusion_matrix_value]
            n += tp
          elsif cell[:confusion_matrix_predicted] == "false" and cell[:confusion_matrix_actual] == "false"
            tn = cell[:confusion_matrix_value]
            n += tn
          elsif cell[:confusion_matrix_predicted] == "false" and cell[:confusion_matrix_actual] == "true"
            fn = cell[:confusion_matrix_value]
            n += fn
          elsif cell[:confusion_matrix_predicted] == "true" and cell[:confusion_matrix_actual] == "false"
            fp = cell[:confusion_matrix_value]
            n += fp
=======
=end
      params[:random_seed] = 2 unless params[:random_seed]
      params[:stratified] = false unless params[:stratified]
      params[:subjectid] = subjectid if subjectid
      uri = OpenTox::RestClientWrapper.post( File.join(CONFIG[:services]["opentox-validation"],"/crossvalidation"),
        params,{:content_type => "text/uri-list"},waiting_task )
      Crossvalidation.new(uri)
    end

    # looks for report for this crossvalidation, creates a report if no report is found
    # @param [String,optional] subjectid
    # @param [OpenTox::Task,optional] waiting_task (can be a OpenTox::Subtask as well), progress is updated accordingly
    # @return [OpenTox::CrossvalidationReport]
    def find_or_create_report( subjectid=nil, waiting_task=nil )
      @report = CrossvalidationReport.find_for_crossvalidation(self, subjectid) unless @report
      @report = CrossvalidationReport.create(self, subjectid, waiting_task) unless @report
      @report
    end

    # PENDING: creates summary as used for ToxCreate
    def summary( subjectid=nil )
      v = YAML.load RestClientWrapper.get(File.join(@uri, 'statistics'),{:accept => "application/x-yaml", :subjectid => subjectid}).to_s
      if v[OT.classificationStatistics]
        res = {
          :nr_predictions => v[OT.numInstances] - v[OT.numUnpredicted],
          :correct_predictions => v[OT.classificationStatistics][OT.percentCorrect],
          :weighted_area_under_roc => v[OT.classificationStatistics][OT.weightedAreaUnderRoc],
        }
        v[OT.classificationStatistics][OT.classValueStatistics].each do |s|
          if s[OT.classValue].to_s=="true"
            res[:true_positives] = s[OT.numTruePositives]
            res[:false_positives] = s[OT.numFalsePositives]
            res[:true_negatives] = s[OT.numTrueNegatives]
            res[:false_negatives] = s[OT.numFalseNegatives]
            res[:sensitivity] = s[OT.truePositiveRate]
            res[:specificity] = s[OT.falsePositiveRate]
            break
#>>>>>>> mguetlein/development:lib/validation.rb
          end
        end
        res
      elsif v[OT.regressionStatistics]
        {
          :nr_predictions => v[OT.numInstances] - v[OT.numUnpredicted],
          :r_square => v[OT.regressionStatistics][OT.rSquare],
          :root_mean_squared_error => v[OT.regressionStatistics][OT.rootMeanSquaredError],
          :mean_absolute_error => v[OT.regressionStatistics][OT.meanAbsoluteError],
        }
      end
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
    # @param [OpenTox::Crossvalidation] 
    # @param [String,optional] subjectid
    # @return [OpenTox::CrossvalidationReport] nil if no report found
    def self.find_for_crossvalidation( crossvalidation, subjectid=nil )
      uris = RestClientWrapper.get(File.join(CONFIG[:services]["opentox-validation"],
        "/report/crossvalidation?crossvalidation="+crossvalidation.uri), {:subjectid => subjectid}).chomp.split("\n")
      uris.size==0 ? nil : CrossvalidationReport.new(uris[-1])
    end
    
    # creates a crossvalidation report via crossvalidation
    # @param [OpenTox::Crossvalidation] 
    # @param [String,optional] subjectid
    # @param [OpenTox::Task,optional] waiting_task (can be a OpenTox::Subtask as well), progress is updated accordingly
    # @return [OpenTox::CrossvalidationReport]
    def self.create( crossvalidation, subjectid=nil, waiting_task=nil )
      uri = RestClientWrapper.post(File.join(CONFIG[:services]["opentox-validation"],"/report/crossvalidation"),
        { :validation_uris => crossvalidation.uri, :subjectid => subjectid }, {}, waiting_task )
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
    # @param [OpenTox::Crossvalidation] 
    # @param [String,optional] subjectid
    # @return [OpenTox::QMRFReport] nil if no report found
    def self.find_for_model( model, subjectid=nil )
      uris = RestClientWrapper.get(File.join(CONFIG[:services]["opentox-validation"],
        "/reach_report/qmrf?model="+model.uri), {:subjectid => subjectid}).chomp.split("\n")
      uris.size==0 ? nil : QMRFReport.new(uris[-1])
    end
    
    # creates a qmrf report via model
    # @param [OpenTox::Model] 
    # @param [String,optional] subjectid
    # @param [OpenTox::Task,optional] waiting_task (can be a OpenTox::Subtask as well), progress is updated accordingly
    # @return [OpenTox::QMRFReport]
    def self.create( model, subjectid=nil, waiting_task=nil )
      uri = RestClientWrapper.post(File.join(CONFIG[:services]["opentox-validation"],"/reach_report/qmrf"), 
        { :model_uri => model.uri, :subjectid => subjectid }, {}, waiting_task )
      QMRFReport.new(uri)
    end
  end
  
end

