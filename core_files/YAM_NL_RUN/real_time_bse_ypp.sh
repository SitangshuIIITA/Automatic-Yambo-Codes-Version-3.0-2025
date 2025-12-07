#!/bin/bash


        cd $DIR
        
        PATH1=$DIR/nscf/$prefix.save
        PATH1=$DIR/Real-Time/RealTime_BSE
        
        
             
                          {
                            echo ""
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Real-time BSE report:"
                           } >> "$DIR/work_report.txt"

        cd $PATH1

                                                        $YPP_NL -u -F ypp_bse.in 2> out
                                                          
							sed -i.bak 's/ETStpsRt= 200/ETStpsRt= 1000/g'                                            "ypp_bse.in"
                                                        sed -i.bak 's/DampMode= "NONE"/DampMode= "LORENTZIAN"/g'                                 "ypp_bse.in"  
                                                        sed -i.bak 's/DampFactor= 0.000000/DampFactor= '"$nl_Damping"'/g'                        "ypp_bse.in" 
                                                        sed -i.bak 's/0.00000 | 20.00000 |/0.00000 | '"$max_NLEnRange"'|/g'                      "ypp_bse.in" 
  
                                                          
                                                          
                                                                            
                                                    	$YPP_NL -F ypp_bse.in 2> out
                                                        
                                                        cd $DIR
                                                        
    
                                            {
                                                        echo "End of real-time BSE calculations."
                                                        echo ""
                                            } >> "$DIR/work_report.txt"


exit;
