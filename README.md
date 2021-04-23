A command line wrapper for the main entry function tapas_physio_main_create_regressors of the Physio (for Philips equipment). 
An additional streamilied version phiwrap.m target mostly Philips equipment.  Parameters and options are to be  provided as command line parameters or in json file. 
User may dump his existing matlab json  physio structure into such file. 
Parameters can be set from  an exisiting example of config, provided by toolset. 
The command line  parameters have highest priority, than ones from the json file.  T
he main purpose of this script is integraton to CBRAIN yet it  can be used with other frameworks as well, or compilation so tool can be used  on machine without MATLAB

 NOTE: All physio-structure can be specified previous to
       running this function, e.g model.retroir.c, 3, save_dir - prefix
       for resulting folder and in_dir are positional parameters, absent 
        in the original physio.

 NEW PARAMETERS
 
 For convenience two positional and two named parameters are added.
 
   in_dir           the name of the folder with input data (no nested)
   
   save_dir         the result will be saved to the folder with name
                    save_dir + in_dir
                    
   presetdefaults   optional preset defaults choice, takes values  
                    'Philips_ECG3T_v1', in the latter 
                    case the parameters will get default values from a corresponding 
                    example script of the toolset;
                    
   paramfile       optional custom file with a json dump of misc parameters, it does
                    overwrite the positional arguments. Note that input or 
                    output files occur in the physio json dump, those will be ignored.
                    either preset defaults or paramfile are recommended
                    otherwise execution will fail
                    
                    also any physio parameters e.g. model.retroivoir.c 
 EXAMPLES

   phiwrap('myfiles_dir', 'results',...
            'paramfile',   'allotherparams.json',...
            'model.retroicor.degree.c', 3) 
