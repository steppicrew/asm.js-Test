
function addHeap( heap, index, data ) {
    var ret= index[0];
    for ( var i= 0; i < data.length; i++ ) {
        heap[index[0]++]= data[i];
    }
    return ret;
}

function run2a( stdlib ) {
    var date= new Date();
    console.log("Start");

    var heap= new ArrayBuffer(128 * 1024);
    // var heap= new Uint8Array(128 * 1024);
    var heap8= heap;
    // var heap8 = new stdlib.Uint8Array(heap);

    var heapIndex= [ 5 ];

    // Tiny Raytracer (C) Gabriel Gambetta 2013
    // ----------------------------------------
    //
    //  Configuration and scene
    //

    var w= 9;

    // Size of the canvas. w is also reused as a "big constant" / "+infinity"
    heap8[0]= addHeap(heap8, heapIndex, [ w ]);

    // Sphere: radius, [cx,  cy,  cz], R,  G,  B, specular exponent, reflectiveness 
    // R, G, B in [0, 9], reflectiveness in [0..9].
    heap8[1] = addHeap(heap8, heapIndex, [
        w,   0, -w, 0,  9, 9, 0,  w,  2,  // Yellow sphere
//        1,   0,  0, 3,  9, 0, 0,  w,  3,  // Red sphere
//        1,  -2,  1, 4,  0, 9, 0,  9,  4,  // Green sphere
//        1,   2,  1, 4,  0, 0, 9,  w,  5,  // Blue sphere
        0
    ]);

    // Ambient light.
    heap8[2] = addHeap(heap8, heapIndex, [ 2 ]);

    // Point lights: intensity, [x,  y,  z]
    // Intensities should add to 10, including ambient.
    heap8[3] = addHeap(heap8, heapIndex, [
        8,  2, 2, 0,
        0
    ]);

    heap8[4] = addHeap(heap8, heapIndex, []);

    var raytracer= Raytracer(window, {}, heap8);
    raytracer.run(0);

    // Get to the raw pixel data.
    var canvas = document.getElementById("c");
    var context2d = canvas.getContext("2d");

    canvas.width = canvas.height = w;

    var image_data= context2d.createImageData(w, w);
    var image_data_data= image_data.data;
    var start= heap8[4];
 console.log(start);
    var i= heap8[start];
    while ( --i >= 0 ) image_data_data[i] = heap8[start + i + 1];

    context2d.putImageData(image_data, 0, 0);

    console.log("End:", new Date() - date);
}

function run2() {
    run2a(window);
}