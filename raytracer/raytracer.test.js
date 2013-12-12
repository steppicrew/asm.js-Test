
function addHeap( heap, index, data ) {
    var ret= index[0];
    for ( var i= 0; i < data.length; i++ ) {
        heap[index[0]++]= data[i];
    }
    return ret;
}

function run2a( stdlib ) {
    var heap= isAsmJs ? new ArrayBuffer(128 * 1024) : [];
    // var heap= new Uint8Array(128 * 1024);
//    var heap8= heap;
    var inHeap= isAsmJs ? new stdlib.Int32Array(heap) : heap;

    var heapIndex= [ 5 ];

    // Tiny Raytracer (C) Gabriel Gambetta 2013
    // ----------------------------------------
    //
    //  Configuration and scene
    //

    var w= 100;

    // Size of the canvas. w is also reused as a "big constant" / "+infinity"
    inHeap[0]= addHeap(inHeap, heapIndex, [ w ]);

    // Sphere: radius, [cx,  cy,  cz], R,  G,  B, specular exponent, reflectiveness 
    // R, G, B in [0, 9], reflectiveness in [0..9].
    inHeap[1] = addHeap(inHeap, heapIndex, [
        w,   0, -w, 0,  9, 9, 0,  w,  2,  // Yellow sphere
        1,   0,  0, 3,  9, 0, 0,  w,  3,  // Red sphere
        1,  -2,  1, 4,  0, 9, 0,  9,  4,  // Green sphere
        1,   2,  1, 4,  0, 0, 9,  w,  5,  // Blue sphere
        0
    ]);

    // Ambient light.
    inHeap[2] = addHeap(inHeap, heapIndex, [ 2 ]);

    // Point lights: intensity, [x,  y,  z]
    // Intensities should add to 10, including ambient.
    inHeap[3] = addHeap(inHeap, heapIndex, [
        8,  2, 2, 0,
        0
    ]);

    inHeap[4] = addHeap(inHeap, heapIndex, []);

    var raytracer= Raytracer(window, {}, heap);
    raytracer.run(0);

    // Get to the raw pixel data.
    var canvas = document.getElementById("c");
    var context2d = canvas.getContext("2d");

    canvas.width = canvas.height = w;

    var image_data= context2d.createImageData(w, w);
    var image_data_data= image_data.data;
    var start= inHeap[4];

    var i= inHeap[start];
console.log('start', start, i)
    var outHeap = isAsmJs ? new stdlib.Uint8Array(heap, (start + 1) * Int32Array.BYTES_PER_ELEMENT ) : heap.slice(start + 1);

console.log('length:', i)
    while ( --i >= 0 ) image_data_data[i] = outHeap[i];

    context2d.putImageData(image_data, 0, 0);

}

function run2() {
    var date= new Date();
    console.log("Start");

    run2a(window);

    console.log("End:", new Date() - date);
}