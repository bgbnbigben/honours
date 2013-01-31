classdef Queue
    %QUEUE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = protected, SetAccess = protected, Hidden = true)
        q
    end
    
    methods
        function [obj] = Queue()
            obj.q = queue_create();
        end
        
        function [out] = size(obj)
           out = queue_size(obj.q);
        end
        
        function push(obj, var)
            queue_push(obj.q, var);
        end
        
        function [out] = pop(obj)
            out = queue_pop(obj.q);
        end
        
        function [out] = peek(obj)
            out = queue_top(obj.q);
        end
        
        function [out] = isempty(obj)
           out = obj.size() == 0;
        end
        
        function delete(obj)
            queue_delete(obj.q);
        end
    end
end

