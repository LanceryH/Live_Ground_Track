let table;
let lon_lat;
let map;
let x;
let y;
let z;

function preload() {
    table = loadTable('satellite.csv', 'csv');
    lon_lat = loadTable('lon_lat.csv', 'csv');
    map = loadImage('MapChart Map.png');
}

function setup() {
    cnv = createCanvas(windowWidth, windowHeight, WEBGL);
    cnv.position(0, 0, "fixed");
    x = table.getRow(0).arr;
    y = table.getRow(1).arr;
    z = table.getRow(2).arr;
}

function draw() {
    // clear();
    background(0);
    orbitControl(1,1,0);
    scale(0.03);

    drawEarth();
    drawTraj();
}

function drawPoint(ind) {
    strokeWeight(10);
    stroke(200,200,0);
    push();
    translate(x[ind]/100, y[ind]/100, z[ind]/100);
    sphere(10);
    pop();
}

function drawTraj() {
    strokeWeight(2);
    stroke(200,0,200);
    push();
    beginShape();
    for (let i = 0; i < x.length; i++) { 
        vertex(x[i]/100, y[i]/100, z[i]/100);
    }
    endShape();
    pop();
}

function drawEarth() {
    noStroke()
    noFill();
    push();
    texture(map)
    rotateX(0);
    rotateY(180);
    rotateZ(0);
    sphere(6e3,50,50);
    pop();
}
