/* global d */
import { select } from 'd3';

import setupData from './data';
import render from './render'; // warning, stats

import addMetricPicker from './toolbar';


export default function init(groupIndex) {
  const data = d[groupIndex];
  const metrics = d.map(d => d.metric);

  // DOM
  const body = select('#main');
  const plotRow = body.insert('div', ':first-child').attr('class', 'viz row');
  const plotDiv = plotRow.append('div').attr('class', 'col-lg-12');
  const controlsRow = plotDiv.append('div').attr('class', 'controls row');
  const svgRow = plotDiv.append('div').attr('class', 'plot row');
  const svgCol = svgRow.append('div').attr('class', 'col-lg-12');
  const svg = svgCol.append('svg');
  const chart = svg.append('g');
  body.insert('h1', ':first-child').text('Alpha Rarefaction');
  chart.append('g').attr('class', 'x axis');
  chart.append('g').attr('class', 'y axis');
  chart.append('text').attr('class', 'x label');
  chart.append('text').attr('class', 'y label');

  // DATA
  const preppedData = setupData(data, data.metric);

  // PLOT
  render(svg, preppedData, data.metric);

  // CONTROLS
  addMetricPicker(controlsRow, metrics, data.metric);

  // STATS
  // stats(body, data);

  // WARNINGS
  // warnings(body, data);
}
