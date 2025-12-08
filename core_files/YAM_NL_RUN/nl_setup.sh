#!/bin/bash


        cd $DIR
        
        PATH1=$DIR/nscf/$prefix.save
        PATH2=$DIR/Real-Time
        
        
             
                          {
                            echo ""
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Nonlinear initialization report:"
                           } >> "$DIR/work_report.txt"

        cd $PATH2

                mkdir -p NL_setup
                mkdir -p $DIR/graph_data/Real_Time_data

        cd $PATH2/NL_setup

                            scp -r $PATH1/SAVE ./.
                                    
                                                        $YAMBO_NL -i -V RL -F setup.in 2> out
                                                        $YAMBO_NL -F setup.in 2> out
                                                        $YPP_NL -y 2> out
                                                                            
 #---------------------------- Define electric field direction vector based on mode_electric_field-----------------------#

                                                        case "$mode_electric_field" in
                                                            x)
                                                                vec="1.000000 | 0.000000 | 0.000000 |"
                                                                ;;
                                                            y)
                                                                vec="0.000000 | 1.000000 | 0.000000 |"
                                                                ;;
                                                            z)
                                                                vec="0.000000 | 0.000000 | 1.000000 |"
                                                                ;;
                                                            xy|yx)
                                                                vec="1.000000 | 1.000000 | 0.000000 |"
                                                               ;;
                                                            yz|zy)
                                                                vec="0.000000 | 1.000000 | 1.000000 |"
                                                               ;;
                                                            zx|xz)
                                                                vec="1.000000 | 0.000000 | 1.000000 |"
                                                                ;;
                                                            *)
                                                                echo "Error: Unsupported mode_electric_field='$mode_electric_field'"
                                                                exit 1
                                                                ;;
                                                        esac

                                                    # Use sed for in-place editing without the backup extension

                                                        sed -i.bak '18s/0.000000 | 0.000000 | 0.000000 |/'"$vec"'/g'                                            "ypp.in"
                                                        sed -i.bak 's/#RmTimeRev/RmTimeRev/g'                                                                   "ypp.in"
                                                        
                                                        $YPP_NL 2> out

                                                        cd FixSymm
                                                        
                                                            $YAMBO_NL -F ../setup.in 2> out
                                                        
                                                        cd $DIR
                                                        
    
                                            {
                                                        echo "End of nonlinear/real-time environment calculations."
                                                        echo ""
                                            } >> "$DIR/work_report.txt"


exit;
