# A wrapper function for testing and timing iterative methods for
# solving the ridge regression problem - 2017 - Cheikh Toure and Robert Gower
function  solveRidge(prob::Prob, method_name, options::MyOptions )

    method = boot_method(method_name,prob,options); 
    if(method=="METHOD DOES NOT EXIST")
        println("FAIL: unknown method name:")
        return
    end
    x = zeros(prob.n)
    d = zeros(prob.n)
    println(method.name);
#     times[1] = toc;
    times= [0];
    if(options.exacterror)
        initial_error= vecnorm(prob.xsol-x);
    end
    initial_residual= vecnorm(prob.A*x -prob.b);
    errors = [1];
    residuals  = [1];
    local tickcounter =1;
    local timeaccum=0;
    iterations =0;
    fail = "failed";
# # Print heard
    if(options.printiters)
         println("-------------------")
         println("It   | Error% | Residual |  Time   ")
         println("-------------------")
    end
    for i = 1:options.maxiter
        tic();
        if(method.name=="CG")
             method.stepmethod(prob,x,d,options,method);
        else
             method.stepmethod(prob,x,d,options);
        end
        x = x+d;
        timeaccum= timeaccum +  toq(); # Keeps track of time accumulated at every iteration
    
        if(mod(i,options.skip_error_calculation)==0 )
             if(options.exacterror)
                    errors= [ errors vecnorm(prob.xsol-x)/initial_error];
             end
             residuals= [residuals vecnorm(prob.A*x -prob.b)/initial_residual];
             times = [ times   timeaccum];
             if(options.printiters)
                ## printing iterations info
                @printf "%3.0d  | %3.2f  |  %3.2f  | %3.4f \n" i 100*errors[end] 100*residuals[end] times[end] ;
             end
             if(errors[end] < options.tol)
                fail ="tol-reached"; iterations =i;
             break;
             end
            if(~isempty(options.restol))
                if(residuals[end] < options.restol)
                    fail ="restol-reached"; iterations =i;
                 break;
            end                
            end
            if(isnan(sum(x)) || isnan(errors[end]) || errors[end] >1000  )  
                 fail = "nan";  iterations = i;
             return;
            end
        end 
     if(timeaccum >options.max_time )
         fail ="times_up";  iterations = i;
         if(options.exacterror)
            errors= [ errors vecnorm(prob.xsol-x)/initial_error];
         end
         residuals= [residuals vecnorm(prob.A*x -prob.b)/initial_residual];
         times = [ times   timeaccum];
         if(options.printiters)
           @printf "%3.0d  | %3.2f  |  %3.2f  | %3.4f \n" i 100*errors[end] 100*residuals[end] times[end] ;
         end
         break;
     end

    end
    if(iterations==0) 
        fail ="max_iter"; 
        iterations=options.maxiter;
    end 
    output = Output(iterations,method.flopsperiter, times, errors, residuals,method.name,fail); 

return output
    
end
