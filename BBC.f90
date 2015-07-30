      PROGRAM BCC_diluido
      implicit none
	  integer, parameter ::L=100, L2=L*L, MCc=200000, m_rand=2147483647
      integer ::i,j,k &
	            passo,iamostras,amostras &
	            Ant(1:L),Suc(1:L) &
                Sigma(-1:1,1:L,1:L,1:L) &
                Bond_i(-1:1,1:L,1:L,1:L) &
                Bond_j(-1:1,1:L,1:L,1:L) &
                Bond_k(-1:1,1:L,1:L,1:L) 
      logical::isolada
      REAL   :: t0,tinc,p &
                W(-8:8,-12:12)& 
                xrand

      OPEN(22,file='hist.dat')
      OPEN(10,access='direct',file='R0072ene.txt')
      OPEN(20,access='direct',file='R0072mag.txt')
      OPEN(30,access='direct',file='R0072hc.txt')
      OPEN(40,access='direct',file='R0072susc.txt')
      OPEN(50,access='direct',file='R0072cumul.txt')
      OPEN(60,access='direct',file='R0072data.txt')
      open(12,file='dados.dat')
     
      
      call Le_Dados(t0,MCc,id,amostra,p)
      MCx = MCc/10
      CALL Inicia_Contorno(L,Ant,Suc)
      CALL Atualiza(t0,W)
       DO iamostras = 1 , amostras
         CALL Inicia_Sigma(L,Sigma,Bond)
         CALL Dilui(L,Bond,m_rand,n_rand,p)
         DO passo = 1 , MCx
            CALL Metropolis(L,L2,Sigma,Bond,W,Ant,Suc,n_rand,m_rand)
         END DO
         DO passo = 1 , MCc
            CALL Metropolis(L,L2,Sigma,Bond,W,Ant,Suc,n_rand,m_rand)
            CALL salvar(passo,L,MCc,Sigma,Bond,Ant,Suc,Dados)
         END DO
         WRITE(*,*) iamostras
      END DO
     
      WRITE(*,*) 'FIM'
      CLOSE(10)
      CLOSE(20)
      CLOSE(30)
      CLOSE(40)
      CLOSE(50)
      CLOSE(60)
 
      END PROGRAM
 
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ok 
      SUBROUTINE Le_Dados(t0,MCc,amostra,p)
      implicit none
      integer MCc, amostra
      double precision t0,p
      read(12,*) p
      read(12,*) t0
      read(12,*) MCc
      read(12,*) amostra
      return
      end
      
      
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ok 
      SUBROUTINE Inicia_Contorno(L,Ant,Suc)
      INTEGER L,i,Ant(1:L),Suc(1:L)
      DO i = 1, L
         Ant(i) = i - 1
         Suc(i) = i + 1
      END DO
      Ant(1) = L
      Suc(L) = 1
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ok
      SUBROUTINE Inicia_Sigma(L,Sigma,Bond)
      INTEGER L,i,j
      INTEGER Sigma(-1:1,1:L,1:L,1:L)
       do sinal =-1,1,2
       DO i = 1 , L
         DO j = 1 , L
           DO k = 1 , L
             Sigma(sinal,i,j,k)=1;
          END DO
         END DO 
        END DO
       END DO 
       RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      SUBROUTINE Dilui(L,Bond_i,Bond_j,Bond_k,id,p)
      INTEGER L,i,j,k,Dilui_count,Dilui_limit, sinal
      INTEGER n_rand,m_rand,aleatorio 
      INTEGER Bond_i(-1:1,1:L,1:L,1:L)
      INTEGER Bond_j(-1:1,1:L,1:L,1:L)
      INTEGER Bond_k(-1:1,1:L,1:L,1:L)
      INTEGER Sigma(-1:1,1:L,1:L,1:L)
      real aux2
      n=0
      DO while(n<2*p*L*L*L)
        i=L*xrand(idum)+1
        j=L*xrand(idum)+1
        k=L*xrand(idum)+1
        aux2=xrand(idum)+0.5
        sinal=aux2/abs(aux2) 
        if (abs(sigma(sinal, i,j,k))==1) then 
          sigma(sinal, i,j,k)=0
	      n=n+1
	    end if 
      end do 
      
      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      SUBROUTINE GERA_LIGACAO(L,Bond_i,Bond_j,Bond_k,id,p)
      INTEGER L,i,j,k,Dilui_count,Dilui_limit, sinal
      INTEGER n_rand,m_rand,aleatorio 
      INTEGER Bond_i(-1:1,1:L,1:L,1:L)
      INTEGER Bond_j(-1:1,1:L,1:L,1:L)
      INTEGER Bond_k(-1:1,1:L,1:L,1:L)
      INTEGER Sigma(-1:1,1:L,1:L,1:L)
      real aux2
      n=0
      sinal=1
       DO i=1,L
         DO j=1,L 
           DO k=1,L 
          if (ISOLADA(sigma(sinal,i,j,k))) then
               Bond_i(-sinal,i,j,k)=1-Bond_i(-sinal,i,j,k)
 	           Bond_i(-sinal,i,sus(j),k)=1 -Bond_i(-sinal,i,sus(j),k)
 	           Bond_i(-sinal,i,j,sus(k))=1 -Bond_i(-sinal,i,j,sus(k))
 	           Bond_i(-sinal,i,sus(j),sus(k))=1 - Bond_i(-sinal,i,sus(j),sus(k))
 	         
 	           Bond_j(-sinal,i,j,k)=1-Bond_j(-sinal,i,j,k)
 	           Bond_j(-sinal,sus(i),j,k)=1- Bond_j(-sinal,sus(i),j,k)
 	     	   Bond_j(-sinal,i,j,sus(k))=1-Bond_j(-sinal,i,j,sus(k))
 	     	   Bond_j(-sinal,sus(i),j,sus(k))=1-Bond_j(-sinal,sus(i),j,sus(k))

  	           Bond_k(-sinal,i,j,k)=1-Bond_k(-sinal,i,j,k)
  	           Bond_k(-sinal,sus(i),j,k)=1-Bond_k(-sinal,sus(i),j,k)
  	           Bond_k(-sinal,i,sus(j),k)=1-Bond_k(-sinal,i,sus(j),k)
  	           Bond_k(-sinal,sus(i),sus(j),k)=1-Bond_k(-sinal,sus(i),sus(j),k)
  	      End if
   end do 
  end do 
 end do   
     sinal=-1
        DO i=1,L
         DO j=1,L 
           DO k=1,L 
            if (ISOLADA(sigma(sinal,i,j,k)) )then
             Bond_i(-sinal,ANT(i),ANT(j),ANT(k))=1 -Bond_i(-sinal,ANT(i),ANT(j),ANT(k))
             Bond_i(-sinal,ANT(i),j,ANT(k))=1 -Bond_i(-sinal,ANT(i),j,ANT(k))
             Bond_i(-sinal,ANT(i),ANT(j),k)=1 -Bond_i(-sinal,ANT(i),ANT(j),k)
 	         Bond_i(-sinal,ANT(i),j,k)=1 - Bond_i(-sinal,i,j,ANT(k))
 	      
 	         Bond_j(-sinal,ant(i),ant(j),ant(k))=1- Bond_j(-sinal,ant(i),ant(j),ant(k))
 	         Bond_j(-sinal,ANT(i),ant(j),k)=1-Bond_j(-sinal,ANT(i),ant(j),k)
 	     	 Bond_j(-sinal,i,ant(j),ANT(k))=1-Bond_j(-sinal,i,ant(j),ANT(k))
 	     	 Bond_j(-sinal,i,ant(j),k)=1-Bond_j(-sinal,i,ant(j),k)

  	         Bond_k(-sinal,ant(i),ant(j),ant(k))=1- Bond_k(-sinal,ant(i),ant(j),ant(k))
  	         Bond_k(-sinal,i,ant(j),ant(k))=1- Bond_k(-sinal,i,ant(j),ant(k))
  	         Bond_k(-sinal,i,j,ant(k))=1-Bond_k(-sinal,i,j,ant(k))
  	         Bond_k(-sinal,ANT(i),j,ant(k))=1-   Bond_k(-sinal,ANT(i),j,ant(k))
  	      End if
        end do 
      end do 
   end do 
      RETURN
      END
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
 FUNCTION  ISOLADA(L,Bond_i,Bond_j,Bond_k,id,p)
	  LOGICAL ISOLADA
      INTEGER L,i,j,k,Dilui_count,Dilui_limit, sinal
      INTEGER n_rand,m_rand,aleatorio
      INTEGER Bond_i(-1:1,1:L,1:L,1:L)
      INTEGER Bond_j(-1:1,1:L,1:L,1:L)
      INTEGER Bond_k(-1:1,1:L,1:L,1:L)
      INTEGER Sigma(-1:1,1:L,1:L,1:L)

      real aux2
        IF(sinal==1)then
           ISOLADA= (ABS(&
                    sigma(-sinal,suc(i),j,k)* &
                    sigma(-sinal,i,suc(j),k)* &
                    sigma(-sinal,suc(i),suc(j),k)* &
                    sigma(-sinal,i,j,suc(k))* &
                    sigma(-sinal,suc(i),j,suc(k))* &
                    sigma(-sinal,i,suc(j),suc(k))* &
                    sigma(-sinal,suc(i),suc(j),suc(k)))==1) .AND. &
                    (sigma(sinal,i,j,k)==0)
         end if
         if(sinal==-1)then
             ISOLADA= (ABS(&
                    sigma(-sinal,i,j,k)* &
                    sigma(-sinal,ant(i),j,k)* &
                    sigma(-sinal,i,ant(j),k)* &
                    sigma(-sinal,i,j,ant(k))* &
                    sigma(-sinal,ant(i),ant(j),k)* &
                    sigma(-sinal,ant(i),j,ant(k))* &
                    sigma(-sinal,i,ant(j),ant(k))* &
                    sigma(-sinal,ant(i),ant(j),ant(k))&
                    )==1) .AND. (sigma(sinal,i,j,k)==0)
           end if                
      RETURN
      END
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE Atualiza(t0,W)
      INTEGER i, j
      REAL t0
      REAL W(-8:8,-12:12) 
      DO i=-8,8
        DO j=-12,12
             AUX=(J1*i+J2*J)/T
             W(i,j)=1
          if (aux >0) W(i,j)=exp(-2*AUX)                    ! probabilidade  METROPOLIS
	    END DO  
	   END DO 
      RETURN
      END
     
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE Salvar(t0,W)
      INTEGER i,j
      INTEGER soma_nj1
      INTEGER soma_nj2
      INTEGER mag
      
	  soma_nj1=0
      soma_nj2=0
	  
	  sinal=1
      DO i=1,L
        DO j=1,L 
          DO k=1,L 
             NJ1= sigma(sinal,i,j,k)*(&
			      sigma(-sinal,i,j,k)+ &
                  sigma(-sinal,suc(i),j,k)+ &
                  sigma(-sinal,i,suc(j),k)+ &
                  sigma(-sinal,suc(i),suc(j),k)+ &
                  sigma(-sinal,i,j,suc(k))+ &
                  sigma(-sinal,suc(i),j,suc(k))+ &
                  sigma(-sinal,i,suc(j),suc(k))+ &
                  sigma(-sinal,suc(i),suc(j),suc(k)) )
                               
	        NJ2= sigma(sinal,i,j,k)*(&
	             sigma(sinal,ant(i),j,k)*Bond_i(sinal,ant(i),j,k)+ &
                 sigma(sinal,suc(i),j,k)*Bond_i(sinal,i,j,k)+ &
                 sigma(sinal,i,ant(j),k)*Bond_j(sinal,i,ant(j),k)+ &
                 sigma(sinal,i,suc(j),k)*Bond_j(sinal,i,j,k)+ &
                 sigma(sinal,i,j,ant(k))*Bond_k(sinal,i,j,ant(k))+ &
                 sigma(sinal,i,j,suc(k))*Bond_k(sinal,i,j,k) )
          
          soma_nj1=NJ1+ soma_nj1
          soma_nj2=NJ2+ soma_nj2
	      END DO 
       END DO  
     END DO  
     
      
     sinal=-1
       DO i=1,L
         DO j=1,L 
           DO k=1,L 
             NJ1= sigma(-sinal,i,j,k)*( &
                  sigma(-sinal,i,j,k)+ &
                  sigma(-sinal,ant(i),j,k)+ & 
                  sigma(-sinal,i,ant(j),k)+ & 
                  sigma(-sinal,ant(i),ant(j),k)+ & 
                  sigma(-sinal,ant(i),j,ant(k))+ &
                  sigma(-sinal,i,j,ant(k))+ &
                  sigma(-sinal,i,ant(j),ant(k))+ &
                  sigma(-sinal,ant(i),ant(j),ant(k)) )  
                               
	        NJ2= sigma(sinal,i,j,k)*  &
	             (sigma(sinal,ant(i),j,k)*Bond_i(sinal,ant(i),j,k)+ &
                 sigma(sinal,suc(i),j,k)*Bond_i(sinal,i,j,k)+ &
                 sigma(sinal,i,ant(j),k)*Bond_j(sinal,i,ant(j),k)+ &
                 sigma(sinal,i,suc(j),k)*Bond_j(sinal,i,j,k)+ &
                 sigma(sinal,i,j,ant(k))*Bond_k(sinal,i,j,ant(k))+ &
                 sigma(sinal,i,j,suc(k))*Bond_k(sinal,i,j,k))
          
           soma_nj1=NJ1+ soma_nj1
           soma_nj2=NJ2+ soma_nj2
	     
	     END DO 
       END DO  
      END DO
     
    DO sinal=-1,1,2
      DO i=1,L
        DO j=1,L 
          DO k=1,L 
             mag= mag+sigma(sinal,i,j,k)
	      END DO 
       END DO  
     END DO
	end do    
     write(20,*) soma_nj1, soma_nj2, mag
    RETURN
    END
     
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE Metropolis(passo,L,MCc,Sigma,Bond,Ant,Suc,Dados)
      INTEGER i, j, k, NJ1, NJ2
      INTEGER Sigma(-1:1,1:L,1:L,1:L)       
	   
	   
	   sinal=1
       DO i=1,L
         DO j=1,L 
           DO k=1,L 
             NJ1= sigma(sinal,i,j,k)*(sigma(-sinal,i,j,k)+ &
                  sigma(-sinal,suc(i),j,k)+ &
                  sigma(-sinal,i,suc(j),k)+ &
                  sigma(-sinal,suc(i),suc(j),k)+ &
                  sigma(-sinal,i,j,suc(k))+ &
                  sigma(-sinal,suc(i),j,suc(k))+ &
                  sigma(-sinal,i,suc(j),suc(k))+ &
                  sigma(-sinal,suc(i),suc(j),suc(k)) )
                               
	        NJ2= sigma(sinal,i,j,k)*(  &
	             sigma(sinal,ant(i),j,k)*Bond_i(sinal,ant(i),j,k)+ &
                 sigma(sinal,suc(i),j,k)*Bond_i(sinal,i,j,k)+ &
                 sigma(sinal,i,ant(j),k)*Bond_j(sinal,i,ant(j),k)+ &
                 sigma(sinal,i,suc(j),k)*Bond_j(sinal,i,j,k)+ &
                 sigma(sinal,i,j,ant(k))*Bond_k(sinal,i,j,ant(k))+ &
                 sigma(sinal,i,j,suc(k))*Bond_k(sinal,i,j,k) )
         
	        IF (XRAND(IDUM)<W(NJ1,NJ2)) then 
			 sigma(sinal,i,j,k)=-sigma(sinal,i,j,k)
			end if      
	     END DO 
       END DO  
     END DO  
     
     
     
     sinal=-1
       DO i=1,L
         DO j=1,L 
           DO k=1,L 
             NJ1= sigma(-sinal,i,j,k)*( &
                  sigma(-sinal,i,j,k)+ &
                  sigma(-sinal,ant(i),j,k)+ & 
                  sigma(-sinal,i,ant(j),k)+ & 
                  sigma(-sinal,ant(i),ant(j),k)+ & 
                  sigma(-sinal,ant(i),j,ant(k))+ &
                  sigma(-sinal,i,j,ant(k))+ &
                  sigma(-sinal,i,ant(j),ant(k))+ &
                  sigma(-sinal,ant(i),ant(j),ant(k)) )  
                               
	        NJ2= sigma(sinal,i,j,k)*(  &
	             sigma(sinal,ant(i),j,k)*Bond_i(sinal,ant(i),j,k)+ &
                 sigma(sinal,suc(i),j,k)*Bond_i(sinal,i,j,k)+ &
                 sigma(sinal,i,ant(j),k)*Bond_j(sinal,i,ant(j),k)+ &
                 sigma(sinal,i,suc(j),k)*Bond_j(sinal,i,j,k)+ &
                 sigma(sinal,i,j,ant(k))*Bond_k(sinal,i,j,ant(k))+ &
                 sigma(sinal,i,j,suc(k))*Bond_k(sinal,i,j,k))
           IF(XRAND(IDUM)<W(NJ1,NJ2)) sigma(sinal,i,j,k)=-sigma(sinal,i,j,k)    
	     END DO  
       END DO  
      END DO
            
      RETURN
      END
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION XRAND(IDUM)
      INTEGER IDUM,IM1,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL XRAND,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211, IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2E-7,RNMX=1.-EPS)
      INTEGER IDUM2,J,K,IV(NTAB),IY
      SAVE IV,IY,IDUM2
      DATA IDUM2/123456719/,IV/NTAB*0/,IY/0/
      IF (IDUM .LE. 0)THEN
      IDUM=MAX(-IDUM,1)
      IDUM2=IDUM
      DO J=NTAB+1,1,-1
      K=IDUM/IQ1
      IDUM=IA1*(IDUM-K*IQ1)-K*IR1
      IF(IDUM .LT. 0)IDUM=IDUM+IM1
      IF(J .LE. NTAB) IV(J)=IDUM
      END DO
      IY=IV(1)
      END IF
      K=IDUM/IQ1
      IDUM=IA1*(IDUM-K*IQ1)-K*IR1
      IF(IDUM .LT. 0)IDUM=IDUM+IM1
      K=IDUM2/IQ2
      IDUM2=IA2*(IDUM2-K*IQ2)-K*IR2
      IF(IDUM2 .LT. 0)IDUM2=IDUM2+IM2
      J= 1 + IY/NDIV
      IY=IV(J)-IDUM2
      IV(J)=IDUM
      IF(IY .LT. 1)IY=IY+IMM1


      XRAND=MIN(AM*IY,RNMX)
      RETURN
      END
