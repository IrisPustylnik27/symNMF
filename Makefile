symnmf: symnmf.o matrix.o
	gcc -ansi -Wall -Wextra -Werror -pedantic-errors symnmf.o matrix.o -o symnmf -lm

symnmf.o: symnmf.c symnmf.h matrix.h
	gcc -ansi -Wall -Wextra -Werror -pedantic-errors -c symnmf.c

matrix.o: matrix.c matrix.h
	gcc -ansi -Wall -Wextra -Werror -pedantic-errors -c matrix.c

clean:
	rm -f *.o symnmf
