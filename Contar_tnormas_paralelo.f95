PROGRAM Contar_tnormas

! Este programa construye y cuenta todas las normas triangulares en una cadena de n+1 elementos
! a partir de un fichero con las t-normas en una cadena de n elementos.

!!!!!!!!!!!!!!!!!
!!! LIBRERÍAS !!!
!!!!!!!!!!!!!!!!!
! Empleamos la librería OpenMP para paralelizar.
USE, INTRINSIC :: omp_lib
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! DECLARACIÓN DE VARIABLES !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL :: start, finish
LOGICAL :: validez
! Para que sólo busque en memoria el espacio para enteros pequeños.
INTEGER, PARAMETER :: tipo = SELECTED_INT_KIND(1)
INTEGER :: dim_ini, num_t_normas_ini, num_t_normas, num_elementos_ini, num_elementos, num_t_normas_final
INTEGER :: i, j, k, p, element, contador_A, IO, dim_vector, dimT
INTEGER :: Tj, Tk, Tjk, mij, Maj
! La variable elemento está escrita en la memoria del TFM como v para simplificar la explicación.
INTEGER(KIND=tipo), ALLOCATABLE:: T(:), A(:,:), elemento(:)

! Para poder leer los ficheros defino las cadenas con los nombres.
CHARACTER(LEN=30) :: filename_matriz, filename_elemento, filename_matriz_nueva, filename_elemento_nuevo, escribir_elemento 
CHARACTER(LEN=30) :: escribir_vector

! Se pide por pantalla la dimensión de la cadena inicial.
PRINT *, 'Dimension de la cadena inicial'
READ(*,*) dim_ini

! Fijamos el tiempo de inicio de la ejecución.
start = OMP_get_wtime()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! LECTURA DE LOS DATOS !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! La longitud del vector A que leemos del caso anterior es la dimensión de la matriz inicial menos el número de elementos de la matriz triangular inferior.
dim_vector = dim_ini-2
DO i = 1, dim_ini-3
    dim_vector = dim_vector + dim_ini -2 - i 
END DO

! En el nuevo vector T tiene que caber el elemento que no es más que la columna nueva que añadimos.
dimT = dim_vector+dim_ini-1
ALLOCATE(T(dimT))
ALLOCATE(elemento(dim_ini-1))

! Escribimos el nombre de los ficheros que vamos a abrir con la dimensión correcta.
! El fichero vectores recoge las t-normas en cadenas de dimensión i0 en forma de vector.
! El fichero elementos guarda los posibles v.
WRITE(filename_matriz ,'("Vectores", i0, ".txt")') dim_ini
WRITE(filename_elemento , '("Elementos_fortran", i0, ".txt")') dim_ini

! Contamos el número de t normas del caso anterior leyendo el fichero.
num_t_normas_ini = 0
OPEN(10, file = filename_matriz, status = 'old')
DO 
    READ(10,*, IOSTAT=IO) 
    IF (IO < 0) EXIT
    num_t_normas_ini = num_t_normas_ini + 1
END DO
CLOSE(10)

! Damos memoria a A que contiene el vector inicial.
ALLOCATE(A(num_t_normas_ini, dim_vector))

! Leemos las matrices y las guardamos:
OPEN(20, file = filename_matriz, status = 'old')
DO i = 1, num_t_normas_ini
    READ(20,*) A(i,:)
END DO
CLOSE(20)

! Contamos los elementos (v en la memoria) que hemos usado en el caso anterior.
num_elementos_ini = 0
OPEN(30, file = filename_elemento, status = 'old')
DO 
    READ(30,*, IOSTAT=IO) 
    IF (IO < 0) EXIT
    num_elementos_ini = num_elementos_ini + 1
END DO
CLOSE(30)


!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! CUERPO DEL PROGRAMA !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Creamos los nombres de los ficheros nuevos.
WRITE(filename_matriz_nueva ,'("Vectores", i0, ".txt")') dim_ini+1
WRITE(filename_elemento_nuevo ,'("Elementos_fortran", i0, ".txt")') dim_ini+1

! Introducimos el formato en el que queremos que se escriba la matriz para ocupar menos espacio en disco.
! Se guarda cada vector en una fila con sus elementos separados por un espacio.
WRITE(escribir_elemento ,'("(",i0, "(I", i1 ,", 1X))")') dim_ini-1, (dim_ini-1)/10 + 1
WRITE(escribir_vector ,'("(",i0, "(I", i1 ,", 1X))")') dim_vector+dim_ini-1, (dim_ini-1)/10 + 1

! Inicializamos los contadores.
num_t_normas = 0
num_t_normas_final = 0

! Abrimos los ficheros que vamos a leer y sobre los que escribiremos en el bucle.
OPEN(40, file = filename_matriz_nueva, status = 'unknown')
OPEN(50, file = filename_elemento_nuevo, status = 'unknown')
OPEN(60, file = filename_elemento, status = 'old')

! Comenzamos la paralelización.
!$OMP PARALLEL PRIVATE(num_t_normas, elemento, p, contador_A, validez, i, j, k, mij, Maj, Tj, Tk, Tjk, T)

! Paralelizamos el bucle que recorre todos los posibles v leyendolos desde el fichero.
!$OMP DO SCHEDULE(STATIC)
DO element = 1, num_elementos_ini
    READ(60,*) elemento(1:dim_ini-2) 
    ! Colocamos estos elementos al final de T.
    T(dim_vector+1:dimT-1) = elemento(1:dim_ini-2) 
    ! Recorremos todos los posibles valores de la última posición de T.
    DO p = elemento(dim_ini-2), dim_ini-1
        elemento(dim_ini-1) = p
        T(dimT) = p 
	! Guardamos el vector elemento una vez le hemos añadido p.
        !$OMP CRITICAL
            WRITE(50, escribir_elemento) elemento
        !$OMP END CRITICAL
	! Recorremos ahora todas las t-normas de la cadena con la dimensión inicial.
        DO contador_A = 1, num_t_normas_ini
	    ! En el momento en que la variable validez sea false, se pasa al siguiente caso.
            validez = .TRUE.
	    ! Este bucle comprueba la monotonía de la t-norma propuesta.
            i = dim_vector - (dim_ini-2) + 1
            DO WHILE (i <= dim_vector)
                IF (A(contador_A, i) > elemento(i - dim_vector + dim_ini-2)) THEN 
                    validez = .FALSE.
                    i = dim_vector
                END IF
                i = i + 1
            END DO
            IF (validez) THEN 
                j = 1
		! A continuación comprobamos la asociatividad.
                ! Separo en 2 el bucle para sólo tener que hacer las comprobaciones en A y elemento y no construir T hasta el final.
                DO WHILE (j <= dim_ini-2)
                    k = 1
                    DO WHILE (k <= j) 
                    ! Si elemento(k) /= 0 entonces elemento(j) /= 0 porque k <= j
                        IF (elemento(k) /= 0)  THEN
                            mij = MIN(elemento(j), k)
                            Maj = MAX(elemento(j), k)
                           ! Vamos a hacer la suma de n números consecutivos con la fórmula n(n+1)/2 y luego sumarle la fila en la que estamos 
                           ! que es siempre el índice más pequeño.
                            Tj = A(contador_A, mij+(Maj-1)*(Maj)/2)
                            Tk = A(contador_A, elemento(k)+(j-1)*j/2)
                            Tjk = A(contador_A, k+(j-1)*j/2)
                            ! El caso en que Tjk es 0 el mismo FORTRAN lo soluciona dando el valor 0 a elemento(Tjk).
                            ! De esta manera no necesitamos tener en cuentra la primera columna de ceros de T.
                            IF ((Tj /= elemento(Tjk)) .OR. (Tk /= elemento(Tjk))) THEN 
                                validez = .FALSE.
                                j = dim_ini
                                k = dim_ini
                            END IF
                        END IF 
                        k = k + 1
                    END DO
                    j = j + 1
                END DO
                
                ! Aqui comprobamos sólamente cuando se recorre la última columna añadida.
                IF (validez) THEN 
                    k = 1
                    DO WHILE (k <= dim_ini-2) 
                        IF (elemento(k) /= 0)  THEN
                            mij = MIN(p, k)
                            Maj = MAX(p, k)
                            IF ((p == dim_ini-1) .AND. (elemento(k) /= elemento(elemento(k)))) THEN
                                validez = .FALSE.
                                k = dim_ini
                            ELSEIF ((p /= dim_ini-1) .AND. (A(contador_A, mij+(Maj-1)*(Maj)/2) /= elemento(elemento(k)))) THEN
                                validez = .FALSE.
                                k = dim_ini
                            END IF
                        END IF 
                        k = k + 1
                    END DO
                END IF 
                IF (validez) THEN 
		    ! Una vez hemos comprobado que la propuesta de t-norma es válida, la guardamos en el fichero
		    ! y sumamos 1 al contador de t-normas.
                    T(1:dim_vector) = A(contador_A, :)
                    num_t_normas = num_t_normas + 1
                    !$OMP CRITICAL
                        WRITE(40, escribir_vector) T
                    !$OMP END CRITICAL 
                END IF 
            END IF
        END DO 
    END DO
END DO
!$OMP END DO NOWAIT
! Para terminar la paralelización, juntamos las t-normas que ha contado cada proceso.
!$OMP CRITICAL
    num_t_normas_final = num_t_normas_final + num_t_normas 
!$OMP END CRITICAL 
!$OMP END PARALLEL
CLOSE(40)
CLOSE(50)
CLOSE(60)

PRINT *, 'Hay ', num_t_normas_final, ' t-normas en la cadena de dimension', dim_ini+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! MOSTRANDO EL TIEMPO !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
finish = OMP_get_wtime()
PRINT '("Tiempo = ",f12.5," segundos.")', finish-start

END PROGRAM Contar_tnormas
