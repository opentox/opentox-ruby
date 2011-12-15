
module OpenTox

  # Class for handling asynchronous tasks
  class Task
    include OpenTox
    attr_accessor :http_code, :due_to_time
    
    def initialize(uri=nil)
      super uri
      @http_code = 202
      @metadata = {
        DC.title => "",
        DC.date => "",
        OT.hasStatus => "Running",
        OT.percentageCompleted => 0.0,
        OT.resultURI => "",
        DC.creator => "", # not mandatory according to API
        DC.description => "", # not mandatory according to API
      }
    end
  
    # Create a new task for the code in the block. Catches halts and exceptions and sets task state to error if necessary. The block has to return the URI of the created resource.
    # @example
    #   task = OpenTox::Task.create do
    #     # this code will be executed as a task
    #     model = OpenTox::Algorithm.run(params) # this can be time consuming
    #     model.uri # Important: return URI of the created resource
    #   end
    #   task.status # returns "Running", because tasks are forked
    # @param [String] title Task title
    # @param [String] creator Task creator
    # @return [OPenTox::Task] Task 
    def self.create( title=nil, creator=nil, max_duration=DEFAULT_TASK_MAX_DURATION, description=nil )
      
      params = {:title=>title, :creator=>creator, :max_duration=>max_duration, :description=>description }
      task_uri = RestClientWrapper.post(CONFIG[:services]["opentox-task"], params, {}, nil, false).to_s
      task = Task.new(task_uri.chomp)

      # measure current memory consumption
=begin
      memory = `free -m|sed -n '2p'`.split
      free_memory = memory[3].to_i + memory[6].to_i # include cache
      if free_memory < 20 # require at least 200 M free memory
        LOGGER.warn "Cannot start task  - not enough memory left (#{free_memory} M free)"
        task.cancel
        return task
        #raise "Insufficient memory to start a new task"
      end

      cpu_load = `cat /proc/loadavg`.split(/\s+/)[0..2].collect{|c| c.to_f}
      nr_cpu_cores = `cat /proc/cpuinfo |grep "cpu cores"|cut -d ":" -f2|tr -d " "`.split("\n").collect{|c| c.to_i}.inject{|sum,n| sum+n}
      nr_cpu_cores = 1 if !nr_cpu_cores
      #if cpu_load[0] > nr_cpu_cores and cpu_load[0] > cpu_load[1] and cpu_load[1] > cpu_load[2] # average CPU load of the last minute is high and CPU load is increasing
      #  LOGGER.warn "Cannot start task  - CPU load too high (#{cpu_load.join(", ")})"
      #  task.cancel
      #  return task
      #  #raise "Server too busy to start a new task"
      #end
=end

      task_pid = Spork.spork(:logger => LOGGER) do
        LOGGER.debug "Task #{task.uri} started #{Time.now}"
        begin
          result = yield task
          LOGGER.debug "Task #{task.uri} done #{Time.now} -> "+result.to_s
          task.completed(result)
        rescue => error
          LOGGER.error "task failed: "+error.class.to_s+": "+error.message
          LOGGER.error ":\n"+error.backtrace.join("\n")
          task.error(OpenTox::ErrorReport.create(error, creator))
        end
      end  
      task.pid = task_pid
      LOGGER.debug "Started task: "+task.uri.to_s
      task
    end  
  
    # Find a task for querying, status changes
    # @param [String] uri Task URI
    # @return [OpenTox::Task] Task object
    def self.find(uri)
      return nil unless uri
      task = Task.new(uri)
      task.load_metadata
      raise "could not load task metadata" if task.metadata==nil or task.metadata.size==0
      task
    end

    # Find a task for querying, status changes
    # @param [String] uri Task URI
    # @return [OpenTox::Task] Task object
    def self.exist?(uri)
      begin
        return find(uri)
      rescue
      end 
    end

    # Get a list of all tasks
    # @param [optional, String] uri URI of task service
    # @return [text/uri-list] Task URIs
    def self.all(uri=CONFIG[:services]["opentox-task"])
      OpenTox.all uri
    end

    def self.from_yaml(yaml)
      @metadata = YAML.load(yaml)
    end
    
    def self.from_rdfxml(rdfxml)
      owl = OpenTox::Parser::Owl.from_rdf(rdfxml, OT.Task)
      task = Task.new(owl.uri)
      task.add_metadata(owl.metadata)
      task
    end

    def to_rdfxml
      s = Serializer::Owl.new
      @metadata[OT.errorReport] = @uri+"/ErrorReport/tmpId" if @error_report
      s.add_task(@uri,@metadata)
      s.add_resource(@uri+"/ErrorReport/tmpId", OT.errorReport, @error_report.rdf_content) if @error_report
      s.to_rdfxml
    end

    def status
      @metadata[OT.hasStatus]
    end

    def result_uri
      @metadata[OT.resultURI]
    end

    def description
      @metadata[DC.description]
    end
    
    def errorReport
      @metadata[OT.errorReport]
    end
    
    def cancel
      RestClientWrapper.put(File.join(@uri,'Cancelled'),{:cannot_be => "empty"})
      load_metadata
    end

    def completed(uri)
      RestClientWrapper.put(File.join(@uri,'Completed'),{:resultURI => uri})
      load_metadata
    end

    def error(error_report)
      raise "no error report" unless error_report.is_a?(OpenTox::ErrorReport)
      RestClientWrapper.put(File.join(@uri,'Error'),{:errorReport => error_report.to_yaml})
      load_metadata
    end
    
    # not stored just for to_rdf
    def add_error_report( error_report )
      raise "not an error report: "+error_report.class.to_s unless error_report.is_a?(ErrorReport)
      @error_report = error_report
    end
    
    def pid=(pid)
      RestClientWrapper.put(File.join(@uri,'pid'), {:pid => pid})
    end

    def running?
      @metadata[OT.hasStatus] == 'Running'
    end

    def queued?
      @metadata[OT.hasStatus] == 'Queued'
    end

    def completed?
      @metadata[OT.hasStatus] == 'Completed'
    end

    def error?
      @metadata[OT.hasStatus] == 'Error'
    end

    def load_metadata
      if (CONFIG[:json_hosts].include?(URI.parse(@uri).host))
        result = RestClientWrapper.get(@uri, {:accept => 'application/x-yaml'}, nil, false)
        @metadata = YAML.load result.to_s
        @http_code = result.code
      else
        @metadata = Parser::Owl::Generic.new(@uri).load_metadata
        @http_code = RestClientWrapper.get(uri, {:accept => 'application/rdf+xml'}, nil, false).code
      end
      raise "could not load task metadata for task "+@uri.to_s if @metadata==nil || @metadata.size==0
    end
    
    # create is private now, use OpenTox::Task.as_task
    #def self.create( params )
      #task_uri = RestClientWrapper.post(CONFIG[:services]["opentox-task"], params, {}, false).to_s
      #Task.find(task_uri.chomp)
    #end
    
=begin
    def self.from_data(data, content_type, code, base_uri)
      task = Task.new(nil)
      task.http_code = code
      task.reload_from_data(data, content_type, base_uri)
      return task
    end

    def reload( accept_header=nil )
      unless accept_header 
        if (CONFIG[:json_hosts].include?(URI.parse(uri).host))
          accept_header = "application/x-yaml"
        else
          accept_header = 'application/rdf+xml'
        end
      end
      result = RestClientWrapper.get(uri, {:accept => accept_header}, false)#'application/x-yaml'})
      @http_code = result.code
      reload_from_data(result, result.content_type, uri)
    end
    
    def reload_from_data( data, content_type, base_uri )
      case content_type
      when /yaml/
        task =  YAML.load data
        TASK_ATTRIBS.each do |a|
          raise "task yaml data invalid, key missing: "+a.to_s unless task.has_key?(a)
          send("#{a.to_s}=".to_sym,task[a])
        end
      when /application\/rdf\+xml/
        owl = OpenTox::Owl.from_data(data,base_uri,"Task")
        self.uri = owl.uri
        (TASK_ATTRIBS-[:uri]).each{|a| self.send("#{a.to_s}=".to_sym, owl.get(a.to_s))}
      else
        raise "content type for tasks not supported: "+content_type.to_s
      end
      raise "uri is null after loading" unless @uri and @uri.to_s.strip.size>0
    end
=end

    # waits for a task, unless time exceeds or state is no longer running
    # @param [optional,OpenTox::Task] waiting_task (can be a OpenTox::Subtask as well), progress is updated accordingly
    # @param [optional,Numeric] dur seconds pausing before cheking again for completion
    def wait_for_completion( waiting_task=nil, dur=0.3)
      
      waiting_task.waiting_for(self.uri) if waiting_task
      due_to_time = Time.new + DEFAULT_TASK_MAX_DURATION
      LOGGER.debug "start waiting for task "+@uri.to_s+" at: "+Time.new.to_s+", waiting at least until "+due_to_time.to_s
      
      load_metadata # for extremely fast tasks
      check_state
      while self.running? or self.queued?
        sleep dur
        #LOGGER.debug "dv ---------------- dur: '#{dur}'"
        dur = dur*2 unless dur>=300.0 
        load_metadata 
        # if another (sub)task is waiting for self, set progress accordingly 
        waiting_task.progress(@metadata[OT.percentageCompleted].to_f) if waiting_task
        check_state
        if (Time.new > due_to_time)
          raise "max wait time exceeded ("+DEFAULT_TASK_MAX_DURATION.to_s+"sec), task: '"+@uri.to_s+"'"
        end
      end
      waiting_task.waiting_for(nil) if waiting_task
      LOGGER.debug "Task '"+@metadata[OT.hasStatus].to_s+"': "+@uri.to_s+", Result: "+@metadata[OT.resultURI].to_s
    end
    
    # updates percentageCompleted value (can only be increased)
    # task has to be running 
    # @param [Numeric] pct value between 0 and 100
    def progress(pct)
      #puts "task := "+pct.to_s
      raise "no numeric >= 0 and <= 100 : '"+pct.to_s+"'" unless pct.is_a?(Numeric) and pct>=0 and pct<=100
      if (pct > @metadata[OT.percentageCompleted] + 0.0001)
        RestClientWrapper.put(File.join(@uri,'Running'),{:percentageCompleted => pct})
        load_metadata
      end
    end
    
    def waiting_for(task_uri)
      RestClientWrapper.put(File.join(@uri,'Running'),{:waiting_for => task_uri})
    end
    
    private
    VALID_TASK_STATES = ["Cancelled", "Completed", "Running", "Error"]
    
    def check_state
      begin
        raise "illegal task state, invalid status: '"+@metadata[OT.hasStatus].to_s+"'" unless 
          @metadata[OT.hasStatus] unless VALID_TASK_STATES.include?(@metadata[OT.hasStatus])
        raise "illegal task state, task is completed, resultURI is no URI: '"+@metadata[OT.resultURI].to_s+
            "'" unless @metadata[OT.resultURI] and @metadata[OT.resultURI].to_s.uri? if completed?
        if @http_code == 202
          raise "#{@uri}: illegal task state, code is 202, but hasStatus is not Running or Queued: '"+@metadata[OT.hasStatus]+"'" unless running? or queued?
        elsif @http_code == 201
          # ignore hasStatus
          # raise "#{@uri}: illegal task state, code is 201, but hasStatus is not Completed: '"+@metadata[OT.hasStatus]+"'" unless completed?
          raise "#{@uri}: illegal task state, code is 201, resultURI is no task-URI: '"+@metadata[OT.resultURI].to_s+
              "'" unless @metadata[OT.resultURI] and @metadata[OT.resultURI].to_s.uri?
        end
      rescue => ex
        raise OpenTox::BadRequestError.new ex.message+" (task-uri:"+@uri+")" 
      end
    end
  end

  # Convenience class to split a (sub)task into subtasks
  #
  # example:
  # a crossvalidation is split into creating datasets and performing the validations
  # creating the dataset is 1/3 of the work, perform the validations is 2/3:
  # Task.as_task do |task|
  #   create_datasets( SubTask.new(task, 0, 33) )
  #   perfom_validations( SubTask.new(task, 33, 100) )
  # end
  # inside the create_datasets / perform_validations you can use subtask.progress(<val>)
  # with vals from 0-100
  #
  # note that you can split a subtask into further subtasks
  class SubTask
    
    def initialize(task, min, max)
      raise "not a task or subtask" if task!=nil and !(task.is_a?(Task) or task.is_a?(SubTask)) 
      raise "invalid max ("+max.to_s+"), min ("+min.to_s+") params" unless 
        min.is_a?(Numeric) and max.is_a?(Numeric) and min >= 0 and max <= 100 and max > min 
      @task = task
      @min = min
      @max = max
      @delta = max - min
    end

    # convenience method to handle null tasks
    def self.create(task, min, max)
      if task
        SubTask.new(task, min, max)
      else
        nil
      end
    end
    
    def waiting_for(task_uri)
      @task.waiting_for(task_uri)
    end
    
    def progress(pct)
      raise "no numeric >= 0 and <= 100 : '"+pct.to_s+"'" unless pct.is_a?(Numeric) and pct>=0 and pct<=100
      #puts "subtask := "+pct.to_s+" -> task := "+(@min + @delta * pct.to_f * 0.01).to_s
      @task.progress( @min + @delta * pct.to_f * 0.01 )
    end
    
    def running?()
      @task.running?
    end
  end

  
  # The David Gallagher feature:
  # a fake sub task to keep the progress bar movin for external jobs
  # note: param could be a subtask
  # 
  # usage (for a call that is normally finished in under 60 seconds):
  #   fsk = FakeSubTask.new(task, 60)
  #   external_lib_call.start
  #   external_lib_call.wait_until_finished
  #   fsk.finished
  #   
  # what happens: 
  # the FakeSubTask updates the task.progress each second until 
  # runtime is up or the finished mehtod is called 
  # 
  # example if the param runtime is too low:
  #   25% .. 50% .. 75% .. 100% .. 100% .. 100% .. 100% .. 100%
  # example if the param runtime is too high:
  #    5% .. 10% .. 15% ..  20% ..  25% ..  30% ..  35% .. 100%
  # the latter example is better (keep the bar movin!) 
  # -> better make a conservative runtime estimate 
  class FakeSubTask
    
    def initialize(task, runtime)
      @task = task
      @thread = Thread.new do
        timeleft = runtime
        while (timeleft > 0 and @task.running?)
          sleep 1
          timeleft -= 1
          @task.progress( (runtime - timeleft) / runtime.to_f * 100 )
        end
      end
    end
    
    # convenience method to handle null tasks
    def self.create(task, runtime)
      if task
        FakeSubTask.new(task, runtime)
      else
        nil
      end
    end
  
    def finished
      @thread.exit
      @task.progress(100) if @task.running?
    end
  end

end
