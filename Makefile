symnmf: symnmf.o matrix.o input.o
	gcc -ansi -Wall -Wextra -Werror -pedantic-errors symnmf.o matrix.o input.o -o symnmf -lm

symnmf.o: symnmf.c symnmf.h matrix.h input.h
	gcc -ansi -Wall -Wextra -Werror -pedantic-errors -c symnmf.c

matrix.o: matrix.c matrix.h
	gcc -ansi -Wall -Wextra -Werror -pedantic-errors -c matrix.c

input.o: input.c input.h matrix.h
	gcc -ansi -Wall -Wextra -Werror -pedantic-errors -c input.c

clean:
	rm -f *.o symnmf
