
// var asm= require('asm.js');
// var validate= asm.validate;
// var ValidationError= asm.ValidationError;

var a= 30;

function fib( n ) {
    if ( n === 0 ) return 0;
    if ( n === 1 ) return 1;
    return fib(n - 1) + fib(n - 2);
};

var start= new Date();
var result= fib(a);
var end= new Date();
document.write('fib(' + a + ')=' + result  + ' in ' + (end-start) + 'ms<br />');

function FibModule( stdlib, foreign, heap ) {
    'use asm';

    function fib2( n ) {
        n= n|0;
        if ( (n|0) == 0 ) return 0;
        if ( (n|0) == 1 ) return 1;
        return ((fib2( ((n|0) - 1)|0 )|0) + (fib2( ((n|0) - 2)|0 )|0))|0;
    }

    return {
        fib: fib2,
    };
};

var Fib= FibModule();

start= new Date();
result= Fib.fib(a);
end= new Date();
document.write('fib_asm(' + a + ')=' + result  + ' in ' + (end-start) + 'ms<br />');

// validate(FibModule.toString());