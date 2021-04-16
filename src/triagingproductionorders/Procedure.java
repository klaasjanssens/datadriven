 /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package triagingproductionorders;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.time.Clock;
import java.util.Random;

public class Procedure {
    /* SET INPUT VALUES */
    int maxPeriod = 200;
    int max_C = 20000; //max number of jobs,customers (can enlarge if necessary)
    int max_run = 10; 
    int max_S = 10; //servers
    int max_AS = 4; //arrival sources (in this case we have 4)
    int max_nr_stations = 10; 
    int max_nr_job_types = 10; 
    int seed;
    
    Random random; // needed to produce random numbers with a given seed
        
    /* COUNTER */
    int i1, i2, i6, run, i3;
    double j1, j2, j3;
    double l1;
    int K, s0, L;
    double [] avg = new double[30];
    char [] naam = new char[300];
    char [] sproblem = new char[10];
    double left_var_Triangular, right_var_Triangular;
    int run_n = 0;
        
    /* INPUT DATA RELATED TO RADIOLOGY DPT */
    int nr_stations;                                                            //Number of workstations
    int [] nr_servers = new int[max_nr_stations];                               //Number of servers per workstation
    
    /* INPUT DATA RELATED TO SYSTEM JOBS */
    int nr_job_types;                                                           // Number of job types
    int nr_workstations_job;                                                    // Number of workstations
    int [] route = new int[max_nr_stations];                                    // Route to follow for each job type
    int [] current_station = new int[max_C];                                    // Matrix that denotes the current station of a job (sequence number)

    
    /* GENERAL DISCRETE EVENT SIMULATION PARAMETERS */
    double t;                                                                   // Simulation time
    int N;                                                                      // Max number of jobs (Stop criterion)
    int T;                                                                      // Max time (Stop criterion)

    /* VARIABLES RELATED TO system JOBS */
    int n;                                                                      // Number of jobs in the system
    int [] n_ws = new int[max_nr_stations];                                     // Number of jobs at a particular workstation
    double [] mean_customers_system = new double[max_run];                      // Throughput time (avg per run)
    double [] tot_n = new double[max_run];                                      // Number of customers in the system over time
    double [][] tot_n_ws = new double[max_run][max_nr_stations];                // Number of customers in a workstation over time
    
    /* PARAMETERS RELATED TO ARRIVAL OF JOBS */
    int nr_arrival_sources;                                                     // Number of arrival sources
    double [] lambda = new double[max_AS];                                      // Arrival rate
    int n_a;                                                                    // Number of jobs arrived to the system
    int [] n_a_ws = new int[max_nr_stations];                                   // Number of jobs arrived to a particular workstation
    double [] t_a = new double[max_AS];                                         // Time of next arrival for each source
    double first_ta;                                                            // First arrival time over all sources;
    int index_arr;                                                              // Source of next arrival;
    double t_lambda;
    double [][] tot_lambda = new double[max_run][max_AS];                       // Counter interarrival time
    int [] job_type = new int[max_C];                                           // Type of job arriving
    double [][] time_arrival = new double[max_run][max_C];                      // Time of arrival of the job to the ancillary services
    double [][][] time_arrival_ws = new double[max_run][max_nr_stations][max_C]; // Time of arrival of a particular job to a particular workstation
    double [] mean_interarrival_time = new double[max_run];
 
    /* PARAMETERS RELATED TO Processing OF JOBS */
    double [][] mu = new double[max_nr_stations][max_nr_job_types];             // Processing rate per station and job type
    double [][] var = new double[max_nr_stations][max_nr_job_types];            // Variance per station and job type
    double [][] sigma = new double[max_nr_stations][max_nr_job_types];          // Standard deviation per station and job type
    int n_d;                                                                    // Number of jobs handled
    int [] n_d_ws = new int[max_nr_stations];                                   // Number of jobs handled in a particular workstation
    int [] n_d_as = new int[max_AS];                                            // Number of jobs processed per type
    double [][] t_d = new double[max_nr_stations][max_S];                       // Time of next departure for each server in each workstation
    double first_td;                                                            // First departure time over all sources
    int index_dep_station;                                                      // Station with first departure
    int index_dep_server;                                                       // Server with first departure
    int index_type;                                                             // Type of next job
    double [] mean_service_time = new double[max_run];                          // Calculated average service time
    double t_mu;                                                                // Generated service time
    double [] tot_mu = new double[max_run];                                     // Total service time generated
    double [][][] time_service = new double[max_run][max_nr_stations][max_C];   // Service time per customer and workstation
    int [][] current_cust = new int[max_nr_stations][max_S];                    // Customer handles by a particular workstation and server
    double [] perc_fail = new double[max_nr_stations];                          // The probability that a job should be re-evaluated at a particular workstation
    int [][] list_process = new int[max_nr_stations][max_C];                    // list of jobs processed at a particular workstation on a particular point in time --> most important!!!
    
    /* PARAMETERS RELATED TO waiting OF JOBS */
    double [] mean_waiting_time = new double[max_run];
    double [] waiting_time = new double[max_run];
    double [][][] waiting_time_job_ws = new double[max_run][max_nr_stations][max_C];    // Waiting time for a job on a particular workstation
    double [] mean_customers_queue = new double[max_run];
    double [] tot_n_queue = new double[max_run];
    double [][] tot_n_queue_ws = new double[max_run][max_nr_stations];          // Total number of jobs in queue at workstation over time
    int [] queue_ws1 = new int[max_C];                           // Total number of jobs in queue at workstation
    int queue_ws1_counter = 0;
    int [] queue_ws2 = new int[max_C];                           // Total number of jobs in queue at workstation
    int queue_ws2_counter = 0;
    int [] queue_ws0 = new int[max_C];                           // Total number of jobs in queue at workstation
    int queue_ws0_counter = 0;
    
    /* VARIABLES RELATED TO Processed JOBS */
    double [][] time_departure = new double[max_run][max_C];
    double [][][] time_departure_ws = new double[max_run][max_nr_stations][max_C];
    double [][] time_system = new double[max_run][max_C];
    double [][][] time_system_job_ws = new double[max_run][max_nr_stations][max_C];
    int [] order_out = new int[max_C];
    double [] mean_system_time = new double[max_run];
    double [][] mean_system_time_as = new double[max_run][max_AS];

    /* OTHER PARAMETERS */
    int infinity;
    double [][][] idle = new double[max_run][max_nr_stations][max_S];
    double [][] rho_ws_s = new double[max_nr_stations][max_S];
    double [] rho_ws = new double[max_nr_stations];
    double rho;
    int [] obj_fct = new int[max_AS];
    double objective;
    int triaging;                                                               // BINARY = 1 if triaging; 0 if no triaging

    /* VARIABLES RELATED TO CLOCK TIME */
    double elapsed_time, time_subproblem; 
    Clock start_time, inter_time, project_start_time;
    
    public Procedure(){
        /* INPUT DATA RELATED TO PRODUCTION DPT */
        nr_stations = 3;
        nr_servers[0] = 2;
        nr_servers[1] = 5;
        nr_servers[2] = 5;

        /* INPUT DATA RELATED TO SYSTEM JOBS */
        nr_job_types = 4;
        triaging  = 0;
        if (triaging == 0){
            nr_workstations_job = 2;
            route[0] = 1;                           
            route[1] = 2;
            route[2] = 0;
        }
        else{
            nr_workstations_job = 3;
            route[0] = 0;                           
            route[1] = 1;
            route[2] = 2;
        }
        
        /* INPUT ARRIVAL PROCESS */
        nr_arrival_sources = 4;
        lambda[0] = 0.5;
        lambda[1] = 1.25;
        lambda[2] = 0.75;
        lambda[3] = 0.2;
        
        /* INPUT SERVICE PROCESS */
        mu[0][0] = 1.0/0.05;
        mu[0][1] = 1.0/0.05;
        mu[0][2] = 1.0/0.05;
        mu[0][3] = 1.0/0.05;
        mu[1][0] = 2.0/3.0;
        mu[1][1] = 2.0/3.0;
        mu[1][2] = 1.0/1.0;
        mu[1][3] = 1.0/1.0;
        mu[2][0] = 4.0/5;
        mu[2][1] = 4.0/5.0;
        mu[2][2] = 2.0/3.0;
        mu[2][3] = 2.0/3.0;
        /**mu[3][0] = 0;
        mu[3][1] = 0;
        mu[3][2] = 0;
        mu[3][3] = 0;*/

        obj_fct[0] = 8;
        obj_fct[1] = 4;
        obj_fct[2] = 2;
        obj_fct[3] = 1;
        
        /* STOP CRITERION (design choice) */
        N = 1000; // Number of jobs
        T = 1000; // Max Time
        
        /* OTHER PARAMETERS */
        infinity = 999999999;
    }
    
    public void doProcedure() throws IOException{
        L = 1;
        
        for (i3 = 0; i3 < L; i3++){
            K = 1;                              //1 replication per run 
            for (run = 0; run < K; run++){
                seed = (i3+1)*K-run;                                            // Ensure you use a different seed each time to get IID replications
                resetVariables();
                run_n = i3;
                production_system();
                output();
            }
        }
        
        /* PRINT OUTPUT of Multiple runs */
        String fileName2 = "Triaging_avg_runs.txt";                             // Output file (CHANGE DIRECTORY)
        File file = new File(fileName2);
        // if file doesnt exists, then create it
	if (!file.exists()) {
            file.createNewFile();                                               // create the file
	} else{
            PrintWriter writer = new PrintWriter(file);                         // empty the file
            writer.print("");
            writer.close();
        }
        FileWriter fileWriter = new FileWriter(file.getAbsoluteFile(),true);    // APPENDS the text file with anything printed to the file during the rest of the procedure
        PrintWriter printWriter = new PrintWriter(fileWriter);                  // OPEN OUTPUT FILE
        
        // TODO STUDENTS

        printWriter.close();
    }
    
    
    private void resetVariables(){ //puts everything to zero 
       /* INPUT DATA RELATED TO SYSTEM JOBS */
        for (i1 = 0; i1 < max_C; i1++)
            current_station[i1] = 0;

        /* GENERAL DISCRETE EVENT SIMULATION PARAMETERS */
        t = 0;
        N = 0;
        
        /* VARIABLES RELATED TO system JOBS */
        n = 0;
        for (i1 = 0; i1 < max_nr_stations; i1++)
            n_ws[i1] = 0;
        
        /* PARAMETERS RELATED TO ARRIVAL OF JOBS */
        n_a = 0;
        first_ta = 0;
        index_arr = 0;
        t_lambda = 0;
        for (i1 = 0; i1 < max_nr_stations; i1++)
            n_a_ws[i1] = 0;
        for (i3 = 0; i3 < max_C; i3++){
            job_type[i3] = 0;
        }
        for (i6 = 0; i6 < max_AS; i6++){
            t_a[i6] = 0;
        }

        /* PARAMETERS RELATED TO Processing OF JOBS */
        for (i1 = 0; i1 < max_nr_stations; i1++){
            n_d_ws[i1] = 0;
            for (i6 = 0; i6 < max_S; i6++){
                t_d[i1][i6] = 0;
                current_cust[i1][i6] = 0;
            }
            for (i6 = 0; i6 < max_C; i6++)
                list_process[i1][i6] = -1;

        }
        n_d = 0;
        first_td = 0;
        index_dep_station = 0;
        index_dep_server = 0;
        t_mu = 0;
        perc_fail[1] = perc_fail[2] = 0;
        
        /* VARIABLES RELATED TO Processed JOBS */
        for (i3 = 0; i3 < max_C; i3++){
            order_out[i3] = 0;
        }

        /* OTHER PARAMETERS */
        rho = 0;
        for (i3 = 0; i3 < max_nr_stations; i3++){
            rho_ws[i3] = 0;
            for (i1 = 0; i1 < max_S; i1++){
                rho_ws_s[i3][i1] = 0;
            }
        }
        
        /* DETERMINE FIRST ARRIVAL + FIRST DEPARTURE */
        
        // TO DO STUDENT    // Put all departure times for all customers to +infty
        
        for (i1 = 0; i1 < nr_stations; i1++){
            for (i2 = 0; i2 < nr_servers[i1]; i2++){   //JUIST?
                t_d[i1][i2] = infinity; 
            }  
        }
      
        // TO DO STUDENT    // Generate first arrival for all sources
        for (i3 = 0; i3 < max_AS; i3++){
            t_a[i3] = Distributions.Poisson_distribution(lambda[i3],this.random); // INVERSION METHOD POISSON DISTRIBUTION
        }
        
        // TO DO STUDENT    // Get next arrival
        for (int i4 = 0; i4 < max_AS; i4++){
            first_ta = infinity;
            if(t_a[i4] < first_ta){
                first_ta = t_a[i4];
                index_arr = i4;
            }
        }
        
        // TO DO STUDENT    // Calculate average arrival time to the system
        
    
    }
    
    private void production_system(){
        
        if(triaging == 0 ){ //Without triaging
            
            // TO DO STUDENT         // Perform simulation until prespecified (time) number of customers have departed (while loop)
             while(n_d < N){                        // As long as the number of customers departed < N, perform simulation
                
                // TO DO STUDENT        //Identify next departure event
                for (i1 = 0; i1 < nr_stations; i1++){
                    for (i2 = 0; i2 < nr_servers[i1]; i2++){
                            first_td = infinity;
                            if(t_d[i1][i2] < first_td){
                             first_td = t_d[i1][i2];
                             index_dep_station = i1;    
                             index_dep_server = i2;
                             //index_type
                            }
                    }
                }
               
                // TO DO STUDENT        // Identify next arrival event
                for (int i4 = 0; i4 < max_AS; i4++){
                    first_ta = infinity;
                    if(t_a[i4] < first_ta){
                    first_ta = t_a[i4];
                    index_arr = i4;
                    }
                }
                
                // TO DO STUDENT        // Identify next event (arrival or departure)
                
                if(first_ta < first_td){ // TO DO STUDENT        // ARRIVAL EVENT
                    arrival_event();
                } else {                // TO DO STUDENT        // DEPARTURE EVENT 
                    departure_event();
                }
                
                
             }
            
        } else { //With triaging
            // TO DO STUDENT        //Identify next departure event
            // TO DO STUDENT        // Identify next arrival event
            // TO DO STUDENT        // Identify next event (arrival or departure)
            // TO DO STUDENT        // ARRIVAL EVENT
            // TO DO STUDENT        // DEPARTURE EVENT
        }
        
        
        
    }
    
    private void arrival_event(){
        // TO DO STUDENT        // DEFINE ARRIVAL EVENT
        
        if(triaging == 0){  //Without triaging
            //WS1 = route[0]
            
            t = first_ta;               // Increment t and jump to arrival time
            
            n++;                        //Increase the total number of jobs in the system
            n_a++;                      //Increase number of jobs arrived to the system
            n_a_ws[1]++;               //Number of jobs arrived to WS1 
            
            job_type[n_a] = index_arr;  // Type of job arriving
            time_arrival[run_n][n_a] = first_ta; //Time of arrival of the nth job for every run
            
           
            int count = 0;
            int worker_idle = 0;
            
            
            for(i1 = 0; i1 < nr_servers[1]; i1++){ //number of workers busy
               if(current_cust[1][i1] ==0){
                count++;
               }              
            }
            
            if(count < nr_servers[1]){
               for(i1 = 0; i1 < nr_servers[1]; i1++){ //first idle worker
                if(current_cust[1][i1] == 0){
                    worker_idle = i1;
                    break;
                }
                                   
            } 
            }
            
            
            if(n_ws[1] < nr_servers[1]){ //Number of jobs at WS1 < number of workers at WS1 --> there is an idle worker/server
                
                time_arrival_ws[run_n][1][n_a] = first_ta; //Initialize arrival time at WS1
                n_ws[1]++;                  //number of jobs at WS1
                
                //Processing of the job
                current_cust[1][worker_idle] = n_a; //Customer handels by WS1 and the idle worker 
                
                
                
                t_mu = Distributions.Exponential_distribution(mu[1][index_arr],this.random);// Generate service time
                time_service[run_n][1][n_a] = t_mu;                                           // Store service time customer n_a
                t_d[1][worker_idle] = t + t_mu;                                             // Generate departure time
                tot_mu[run_n] += t_mu;                                                        //  Update Total Service Time
                
                current_station[n_a] = route[0];                                                   //Current station of a job 
                
            } else {                        //In queue
                //Add unit to queue
                int c = 0;
                while(c == 0){
                    for(int q = 0; q <queue_ws1.length;q++ ){
                        if(queue_ws1[q] ==0){
                            queue_ws1[q] = n_a;
                            c++;
                        }
                    }
                }
                queue_ws1_counter++;
                
            }
            
            // Generate interarrival time of next arrival
            for (i3 = 0; i3 < max_AS; i3++){
            t_a[i3] = Distributions.Exponential_distribution(lambda[i3],this.random); 
            } 
        
            // Calculate arrival time of next arrivals (t+ta)
            for (i3 = 0; i3 < max_AS; i3++){
            t_a[i3] = t_a[i3] + t;
            }
                
            // mean_interarrival_time[i3]
            
        } else {            //With triaging
            
        }
        
    }

    private void departure_event(){
        // TO DO STUDENT        // DEFINE DEPARTURE EVENT
        
        //In which station is the departure event?
        //WS1, WS2 or WS3?
        
        
        //WS1                     
        if(index_dep_station == 1){
           
            
            n_d_ws[index_dep_station]++;        //number of jobs handels at WS1
            t = first_td;       //update time
            
            //1. Idle maken unit eruit
            n_ws[index_dep_station]--;                  //number of jobs at WS1
            
                            //List of jobs processed at a particular WS on a particular moment in time - niet zeker aan wat gelijk
            list_process[1][n_d_ws[1]] = current_cust[index_dep_station][index_dep_server]; //n_a opslaan            
            current_cust[1][index_dep_server] = 0; 
            
            //2. Processed 
            time_departure_ws[run_n][index_dep_station][list_process[1][n_d_ws[1]]] = t; //Save departure time ws
                //Time in WS1 
            time_system_job_ws[run_n][index_dep_station][list_process[1][n_d_ws[1]]] = t - time_arrival_ws[run_n][index_dep_station][list_process[1][n_d_ws[1]]];
                //Waiting time queue WS1
            waiting_time_job_ws[run_n][index_dep_station][list_process[1][n_d_ws[1]]] = time_arrival_ws[run_n][index_dep_station][list_process[1][n_d_ws[1]]] - time_arrival[run_n][n_d_ws[1]];
            
            //3. Units queue?
            
            if( queue_ws1_counter > 0){ //Start processing of next unit in WS1
                int next_arrival = 0;
                
                queue_ws1_counter--;
                next_arrival = queue_ws1[0];
                
                
                for(int q = 1; q <queue_ws1.length;q++ ){
                    queue_ws1[q-1] = queue_ws1[q];
                    }
                
            
                //Last post fill again
                queue_ws1[max_C - 1] =0;
                
                //Processing in WS1
                time_arrival_ws[run_n][index_dep_station][next_arrival] = first_ta; //Initialize arrival time at WS1
                n_ws[1]++;                  //number of jobs at WS1
                
                //Processing of the job
                current_cust[index_dep_station][index_dep_server] = next_arrival; //Customer handels by WS1 and the idle worker 
                
                t_mu = Distributions.Exponential_distribution(mu[index_dep_station][job_type[next_arrival]],this.random);// Generate service time
                time_service[run_n][index_dep_station][next_arrival] = t_mu;                                            // Store service time customer n_a
                t_d[index_dep_station][index_dep_server] = t + t_mu;                                                    // Generate departure time
                tot_mu[run_n] += t_mu;                                                                                  //  Update Total Service Time
                
                current_station[next_arrival] = route[0];                                                               //Current station of a job 
            }   
            } else {
                t_d[index_dep_station][index_dep_server] = infinity;
            }
            
            
        }
        //WS2 - unit is processed
        else if(index_dep_station == 2){
            int WS2 = 2;
            
            
            //1. Update statistics
            n--;                        //Decrease the total number of jobs in the system
            //tot_n[i3]--;               //Decrease number of customers in the system over time
            
            n_d++;                       //Update number of units departed from the system
            t = first_td;                                            // Increment t and jump to departure time
            
            //2. Update waiting queue
            if (n_a_ws[WS2] >= 1){                                    // Update total number of customers in queue of WS2
                n_a_ws[WS2]--;                                        //Remove unit from the queue --> unit will be processed in W 
            }
            
            //time_system += (t-time_arrival[n_d]);               // Update time in system (overall)
            time_departure_ws[i3][WS2][n_d] = t;                   // Store departure time of customer n_d in WS2
            time_departure[i3][n_d] = t;                           // Store deapurture time of customer n_d 
             
             
            if (n_a_ws[WS2] == 0){                                        // Set departure time +inf if no customers are left in the queue 
                   t_d[WS2][index_dep_server] =infinity;                    // only the server!
                } else{
                    t_mu = Distributions.Exponential_distribution(mu,this.random);// If there are still customers in the queue, generate service time for next customer
                    time_service[n_d+1] = t_mu;                     // Store service for next customer
                    t_d = t + t_mu; tot_mu += t_mu;                 // Change next departure time
                }
            
            

            
        }
        //WS0 -triaging
        else if(index_dep_station == 0){
            
        }
        
        
        
    }
    
    private void output() throws IOException{
        String fileName1 = "Output_Triaging" + run + ".txt";
        File file = new File(fileName1);
        // if file doesnt exists, then create it
        if (!file.exists()) {
            file.createNewFile();                                               // create the file
        } else {
            PrintWriter writer = new PrintWriter(file);                         // empty the file
            writer.print("");
            writer.close();
        }
        FileWriter fileWriter = new FileWriter(file.getAbsoluteFile(), true);   // APPENDS the text file with anything printed to the file during the rest of the procedure
        PrintWriter printWriter = new PrintWriter(fileWriter);                  // OPEN OUTPUT FILE
        
        for (i1=0; i1<nr_stations; i1++){                                       // PRINT Utilisation
            printWriter.println("Utilisation servers Station WS" + i1 + ":\t");
            for(i2=0;i2<nr_servers[i1]; i2++){
                j1 = (idle[run][i1][i2]/t);
                rho_ws_s[i1][i2] = 1-j1;
                printWriter.println(rho_ws_s[i1][i2]);
            }
            printWriter.println("\n");
        }
        printWriter.println("\n");

        for (i1=0; i1<nr_stations; i1++){
            printWriter.println("Avg utilisation Station WS" + i1 + ":\t");
            for(i2=0;i2<nr_servers[i1]; i2++){
                rho_ws[i1] += rho_ws_s[i1][i2];
            }
            rho_ws[i1] = rho_ws[i1]/nr_servers[i1];
            printWriter.println(rho_ws[i1] + "\n");
        }
        printWriter.println("\n");

        for (i1=0; i1<nr_stations; i1++){
           rho += rho_ws[i1];
        }
        rho /= nr_stations;
        printWriter.println("Overall avg utilisation:" + rho + "\n\n");

        printWriter.println("\n");

        for (i1 = 0; i1 < n_d; i1++){                                           // PRINT system time = cycle time (observations and running average)
            mean_system_time[run] += time_system[run][order_out[i1]];
        }
        printWriter.println("Cycle time\n\n");
        j1 = mean_system_time[run]/n_d;
        printWriter.println("Avg cycle time:" + j1 + "\n\n");

        mean_system_time[run] = 0;
        for (i1 = 0; i1 < nr_arrival_sources; i1++){
            mean_system_time_as[run][i1] = 0;
            n_d_as[i1] = 0;
        }
        printWriter.println("Number\tObservation\tRunning Average\n");

        for (i1 = 0; i1 < n_d; i1++){                                           // Calculate cycle time per order type
            mean_system_time[run] += time_system[run][order_out[i1]];
            mean_system_time_as[run][job_type[order_out[i1]]] += time_system[run][order_out[i1]];
            n_d_as[job_type[order_out[i1]]]++;
            j1 = mean_system_time[run]/(i1+1);
            printWriter.println(i1 + "\t" + time_system[run][order_out[i1]] + "\t" + j1);
        }

        printWriter.println("Arr S\tNumber\tCycle time\n");
        objective = 0;
        for (i1 = 0; i1 < nr_arrival_sources; i1++){                            // PRINT cycle time per order type
            if (n_d_as[i1]==1){
                j1 = mean_system_time_as[run][i1]/n_d_as[i1];
            }
            else{
                j1 = 0;
            }
            printWriter.println(i1 + "\t" + n_d_as[i1] + "\t" + j1);
            objective += obj_fct[i1] * j1;                                      //mean_system_time_as[run][i1];
        }

        printWriter.println("Objective: " + objective + "\n");
        printWriter.close();
    }
}
