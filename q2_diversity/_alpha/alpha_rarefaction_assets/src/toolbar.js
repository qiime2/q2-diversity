/* global d */
/* global document */
/* global XMLSerializer */

// import { select } from 'd3';

// import setupData from './data';
// import { render, stats, warnings } from './render';


export function addMetricPicker(row, metrics, selectedMetric) {
  const grp = row.append('div').attr('class', 'col-lg-2 form-group metricPicker');
  grp.append('label').text('Metric');
  grp.append('select')
    .attr('class', 'form-control')
    // .on('change', function changeCategory() {
    //   const data = d[this.selectedIndex];
    //   const svg = select('svg');
    //   const preppedData = setupData(data);
    //   render(svg, preppedData);
    //   const body = select('body .container-fluid');
    //   stats(body, data);
    //   warnings(body, data);
    // })
    .selectAll('option')
    .data(metrics)
    .enter()
      .append('option')
      .attr('value', d => d)
      .text(d => d)
      .property('selected', d => (d === selectedMetric));
  return grp;
}

export function addDownloadLinks(sel, metric) {
  const grp = sel.append('div').attr('class', 'col-lg-2 form-group');
  grp.append('label').html('&nbsp;');
  grp.append('button')
    .text('Download CSV')
    .attr('class', 'btn btn-default form-control')
    .on('click', () => {
      const url = `metric-${metric}.csv`;
      const link = document.createElement('a');
      link.setAttribute('href', url);
      link.setAttribute('download', url);
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
    });
}
