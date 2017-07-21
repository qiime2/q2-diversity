import {
  scaleLinear,
  axisBottom,
  axisLeft,
  // line,
  // curveCardinal,
  // scaleOrdinal,
  // keys,
  // schemeCategory10,
} from 'd3';

import { setupXLabel, setupYLabel } from './axis';

// function renderPlot(svg, data, x, y) {
  // const chart = svg.select('g');
  // chart.selectAll('path').remove();
  // chart.selectAll('dot').remove();

  // const depthIndex = data.data.columns.indexOf('depth');
  // const medianIndex = data.data.columns.indexOf('median');
  // // const sampleIdIndex = data.data.columns.indexOf('sample-id');

  // const valueline = line()
  //   .x(d => x(d[depthIndex]))
  //   .y(d => y(d[medianIndex]))
  //   .curve(curveCardinal);

  // const points = [data.data.data.sort((a, b) => a[depthIndex] - b[depthIndex])];

  // const color = scaleOrdinal(schemeCategory10)
  //   .domain(keys(data[sampleIdIndex]).filter());

  // chart.append('path')
  //   .data(points)
  //   .attr('class', 'line')
  //   .style('fill', 'none')
  //   .style('stroke', 'blue')
  //   .attr('d', valueline);

  // chart.selectAll('dot')
  //     .data(points)
  //   .enter()
  //     .append('circle')
  //       .attr('cx', d => x(d[depthIndex]))
  //       .attr('cy', d => y(d[medianIndex]))
  //       .attr('r', 4)
  //       .style('stroke', 'green');

  // const samples = chart.selectAll('.sample')
  //     .data(points)
  //   .enter().append('g')
  //     .attr('class', 'sample');

  // samples.append('path')
  //     .attr('class', 'line')
  //     .attr('d', d => valueline(d.values));
      // .style('stroke', d => color(d[sampleIdIndex]));
// }

export default function render(svg, data) {
  const height = 400;
  const width = 1000;
  const margin = { top: 20, left: 70, right: 50, bottom: 50 };
  const chart = svg.select('g');

  const { xAxisLabel, yAxisLabel, minX, maxX, minY, maxY } = data;

  const xAxis = axisBottom();
  const yAxis = axisLeft();

  let pad = (maxX - minX) * 0.03;
  if (Number.isInteger(minX) && Number.isInteger(maxX)) {
    pad = Math.max(Math.round(pad), 1);
    const between = Math.max(3, (maxX - minX) + (2 * pad));
    xAxis.ticks(Math.min(between, 12), 'd');
  }

  const x = scaleLinear().domain([minX - pad, maxX + pad]).range([0, width]).nice();
  const y = scaleLinear().domain([minY, maxY]).range([height, 0]).nice();

  xAxis.scale(x);
  yAxis.scale(y);

  chart.attr('transform', `translate(${margin.left},${margin.top})`);

  setupXLabel(svg, width, height, xAxisLabel, xAxis);
  setupYLabel(svg, height, yAxisLabel, yAxis);

  // renderPlot(svg, data, x, y);

  svg.attr('width', width + margin.left + margin.right)
    .attr('height', height + margin.bottom + margin.top);
}

// export function stats(body, data) {
//   const { stats: { method, testStat, pVal, sampleSize } } = data;
//   select('#method').text(method);
//   select('#test-stat').text(testStat);
//   select('#p-val').html(pVal);
//   select('#sample-size').html(sampleSize);
// }
