  r  D   k820309              10.1        ÿóL                                                                                                           
       splines.f90 SPLINES              SPLINE_TYPE SPLINE_NULL SPLINE_INIT SPLINE_END SPLINE_EVAL SPLINE_EVAL_DERIV SPLINE_EVAL_DERIV2 SPLINE_EVAL_INTEG SPLINE_INTEG SPLINE_DERIV SPLINE_DERIV2 SPLINE_MESH_TRANSFER i@| i@                                            
                                                  
                                                         |  #COPY_SPLINES    #         @     @X                                       #SPLINE_A    #SPLINE_B              
D @@                                          #SPLINE_TYPE              
   @                                         #SPLINE_TYPE                                                       #EQUAL_SPLINES    %         @   @X                                             #EQUAL_SPLINES%ASSOCIATED    #SPLINE_A 	   #SPLINE_B 
                 @                         ASSOCIATED           
   @                      	                   #SPLINE_TYPE              
   @                      
                   #SPLINE_TYPE                      @                        '                    #INFO                 D                                                 #SPLINE_INFO_TYPE                  @  @                       '                    #ACC    #SPL    #N_PTRS                 D                                                     D                                                    D                                        #         @                                            #SPLINE              D  @                                          #SPLINE_TYPE    #         @                                            #SPLINE    #N    #X    #Y    #INTERP_TYPE              
D @@                                          #SPLINE_TYPE              
  @@                                          
   @                                         
    p          5  p        r        5  p        r                               
   @                                         
    p          5  p        r        5  p        r                                
  @@                                 #         @                                           #SPLINE_END%ASSOCIATED    #SPLINE                  @                         ASSOCIATED           
D @@                                          #SPLINE_TYPE    %         @                                         
       #SPLINE    #X              
   @                                         #SPLINE_TYPE              
  @@                          
      %         @                                         
       #SPLINE    #X               
   @                                         #SPLINE_TYPE              
  @@                           
      %         @                      !                   
       #SPLINE "   #X #             
   @                      "                   #SPLINE_TYPE              
  @@                     #     
      %         @                      $                   
       #SPLINE %   #A &   #B '             
   @                      %                   #SPLINE_TYPE              
  @@                     &     
                
  @@                     '     
      %         @                       (                   
       #NP )   #X *   #Y +   #A ,   #B -   #INTERP_TYPE .             
  @@                      )                    
  @@                     *                    
    p          5  p        r )       5  p        r )                              
  @@                     +                    
    p          5  p        r )       5  p        r )                               
  @@                     ,     
                
  @@                     -     
                
  @@                      .           (        `                       /                   
                
    #NP 0   #X 1   #Y 2   #INTERP_TYPE 3   p          5 O p            5 O p                                    
  @@                      0                    
  @@                     1                    
    p          5  p        r 0       5  p        r 0                              
  @@                     2                    
    p          5  p        r 0       5  p        r 0                               
  @@                      3           (        `                       4                                   
    #NP 5   #X 6   #Y 7   #INTERP_TYPE 8   p          5 O p            5 O p                                    
  @@                      5                    
  @@                     6                    
    p          5  p        r 5       5  p        r 5                              
  @@                     7                    
    p          5  p        r 5       5  p        r 5                               
  @@                      8           #         @                          9                   #NP_A :   #X_A ;   #Y_A <   #NP_B =   #X_B >   #Y_B ?   #INTERP_TYPE @             
  @@                      :                    
  @@                     ;                    
    p          5  p        r :       5  p        r :                              
  @@                     <                    
    p          5  p        r :       5  p        r :                               
   @                      =                    
  @@                     >                    
    p          5  p        r =       5  p        r =                              D  @                     ?                    
     p          5  p        r =       5  p        r =                               
  @@                      @                        fn#fn    ¼   È   b   uapp(SPLINES      8   J  GLOBAL    ¼  8   J  GSL_INTERFACE    ô  J      i@|    >  \      COPY_SPLINES &     Q   a   COPY_SPLINES%SPLINE_A &   ë  Q   a   COPY_SPLINES%SPLINE_B    <  K      i@            EQUAL_SPLINES )   	  ;      EQUAL_SPLINES%ASSOCIATED '   D  Q   a   EQUAL_SPLINES%SPLINE_A '     Q   a   EQUAL_SPLINES%SPLINE_B    æ  R       SPLINE_TYPE !   8  ^   !   SPLINE_TYPE%INFO !     f       SPLINE_INFO_TYPE %   ü  @   !   SPLINE_INFO_TYPE%ACC %   <  @   !   SPLINE_INFO_TYPE%SPL (   |  @   !   SPLINE_INFO_TYPE%N_PTRS    ¼  L       SPLINE_NULL #     Q   a   SPLINE_NULL%SPLINE    Y  r       SPLINE_INIT #   Ë  Q   a   SPLINE_INIT%SPLINE      8   a   SPLINE_INIT%N    T  ¬   a   SPLINE_INIT%X     	  ¬   a   SPLINE_INIT%Y (   ¬	  8   a   SPLINE_INIT%INTERP_TYPE    ä	  g       SPLINE_END &   K
  ;      SPLINE_END%ASSOCIATED "   
  Q   a   SPLINE_END%SPLINE    ×
  [       SPLINE_EVAL #   2  Q   a   SPLINE_EVAL%SPLINE      8   a   SPLINE_EVAL%X "   »  [       SPLINE_EVAL_DERIV )     Q   a   SPLINE_EVAL_DERIV%SPLINE $   g  8   a   SPLINE_EVAL_DERIV%X #     [       SPLINE_EVAL_DERIV2 *   ú  Q   a   SPLINE_EVAL_DERIV2%SPLINE %   K  8   a   SPLINE_EVAL_DERIV2%X "     b       SPLINE_EVAL_INTEG )   å  Q   a   SPLINE_EVAL_INTEG%SPLINE $   6  8   a   SPLINE_EVAL_INTEG%A $   n  8   a   SPLINE_EVAL_INTEG%B    ¦  }       SPLINE_INTEG     #  8   a   SPLINE_INTEG%NP    [  ¬   a   SPLINE_INTEG%X      ¬   a   SPLINE_INTEG%Y    ³  8   a   SPLINE_INTEG%A    ë  8   a   SPLINE_INTEG%B )   #  8   a   SPLINE_INTEG%INTERP_TYPE    [  Ó       SPLINE_DERIV     .  8   a   SPLINE_DERIV%NP    f  ¬   a   SPLINE_DERIV%X      ¬   a   SPLINE_DERIV%Y )   ¾  8   a   SPLINE_DERIV%INTERP_TYPE    ö  Ó       SPLINE_DERIV2 !   É  8   a   SPLINE_DERIV2%NP       ¬   a   SPLINE_DERIV2%X     ­  ¬   a   SPLINE_DERIV2%Y *   Y  8   a   SPLINE_DERIV2%INTERP_TYPE %            SPLINE_MESH_TRANSFER *     8   a   SPLINE_MESH_TRANSFER%NP_A )   R  ¬   a   SPLINE_MESH_TRANSFER%X_A )   þ  ¬   a   SPLINE_MESH_TRANSFER%Y_A *   ª  8   a   SPLINE_MESH_TRANSFER%NP_B )   â  ¬   a   SPLINE_MESH_TRANSFER%X_B )     ¬   a   SPLINE_MESH_TRANSFER%Y_B 1   :  8   a   SPLINE_MESH_TRANSFER%INTERP_TYPE 