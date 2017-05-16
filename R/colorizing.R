colorizing <- function(excm){
  color="black"
    if(excm > .75)
      {color="mediumpurple4"
      }else{if(excm > .5)
        {color="green3"
        }else{if(excm > .25)
          {color="turquoise2"
          }else{if(excm > .05)
            {color="blue"
            }
          }
        }
      }   
  return(color)
}